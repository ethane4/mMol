// MovieExportSheet.swift — render the Timeline to a movie file with in-app
// encoding (H.264 / HEVC / ProRes via AVAssetWriter, GIF via ImageIO, or a PNG
// sequence) — no ffmpeg needed.
//
// Pipeline: pause playback, apply the export-only quality overrides (#581),
// then per frame set cmd.frame(N) and capture the full Metal pipeline offscreen
// via engine.renderHiResPNG (the known-good capture path — NOT cmd.mpng, which
// uses a GL framebuffer path unavailable on this backend), load the PNG,
// downsample it if supersampling, and stream it into the encoder. Frames are
// written/released one at a time so memory stays bounded. The overrides are
// restored when the export completes, is cancelled or fails. The result is
// saved through a Save panel.

import AVFoundation
import AppKit
import CoreGraphics
import ImageIO
import SwiftUI
import UniformTypeIdentifiers

// MARK: - Options

enum MovieQuality: String, Codable, CaseIterable, Identifiable {
  case draft, standard, high, maximum, custom
  var id: String { rawValue }
  var label: String { rawValue.prefix(1).uppercased() + rawValue.dropFirst() }

  // Export-only setting overrides per preset. KEEP IN SYNC with
  // QUALITY_PRESETS in modules/pymol/movie_exporting.py (cmd.movie_export),
  // so scripted and UI exports at the same preset match. `standard` is empty
  // on purpose: no overrides = exactly the pre-#581 export.
  var overrides: [String: Int] {
    switch self {
    case .draft:
      return ["metal_msaa": 0]
    case .standard, .custom:
      return [:]
    case .high:
      return [
        "cartoon_sampling": 14, "cartoon_oval_quality": 16,
        "cartoon_tube_quality": 12, "cartoon_loop_quality": 10,
        "ribbon_sampling": 4, "sphere_quality": 3, "stick_quality": 16,
        "cgo_sphere_quality": 3, "metal_msaa": 1,
      ]
    case .maximum:
      return [
        "cartoon_sampling": 20, "cartoon_oval_quality": 24,
        "cartoon_tube_quality": 16, "cartoon_loop_quality": 12,
        "ribbon_sampling": 8, "sphere_quality": 4, "stick_quality": 24,
        "cgo_sphere_quality": 4, "surface_quality": 2,
        "metal_msaa": 1,
      ]
    }
  }

  // Mirrors PRESET_SUPERSAMPLE in movie_exporting.py.
  var supersample: Int { self == .maximum ? 2 : 1 }
}

struct MovieExportOptions: Codable, Equatable {
  enum Format: String, Codable, CaseIterable, Identifiable {
    case video = "Video"
    case gif = "GIF"
    case png = "PNG sequence"
    var id: String { rawValue }
  }
  enum Codec: String, Codable, CaseIterable, Identifiable {
    case h264 = "H.264"
    case hevc = "HEVC"
    case prores = "ProRes 422"
    var id: String { rawValue }
  }

  var format: Format = .video
  var codec: Codec = .h264
  var width = 1280
  var height = 720
  var quality: MovieQuality = .standard
  var overrides: [String: Int] = [:]  // what's applied; preset's table unless .custom
  var supersample = 1  // 1, 2 or 4: render larger, downsample
  var fpsOverride = 0  // 0 = the timeline's movie_fps
  var bitrateMbps = 0.0  // 0 = let AVFoundation choose (H.264/HEVC)
  var movContainer = false  // .mov even for H.264/HEVC (ProRes always .mov)

  var fileExtension: String {
    switch format {
    case .gif: return "gif"
    case .png: return ""
    case .video: return (codec == .prores || movContainer) ? "mov" : "mp4"
    }
  }

  // Largest offscreen render the exporter will attempt (render = output ×
  // supersample). Bigger targets than this risk the GPU's texture limit and
  // gigabytes of MSAA intermediates.
  static let maxRenderDimension = 8192

  func allows(supersample s: Int) -> Bool {
    max(width, height) * s <= Self.maxRenderDimension
  }

  mutating func applyPreset(_ q: MovieQuality) {
    quality = q
    guard q != .custom else { return }
    overrides = q.overrides
    // Keep the preset's factor even when it exceeds the render limit at this
    // size: the exporter falls back to 1× and the sheet warns about it.
    supersample = q.supersample
  }
}

// A cmd.movie_export request (the JSON movie_exporting.movie_export prints).
struct MovieExportRequest {
  var path: String
  var options: MovieExportOptions
  var first: Int
  var last: Int
  var frames: Int  // count_frames() when the request was made

  private struct Wire: Decodable {
    var path: String
    var frames: Int
    var width: Int
    var height: Int
    var format: String
    var codec: String
    var quality: String
    var overrides: [String: Int]
    var supersample: Int
    var fps: Int
    var bitrate: Double
    var first: Int
    var last: Int
  }

  struct DecodeError: Error { let message: String }

  static func decode(_ json: String) -> Result<MovieExportRequest, DecodeError> {
    func bad(_ m: String) -> Result<MovieExportRequest, DecodeError> {
      .failure(DecodeError(message: m))
    }
    guard let w = try? JSONDecoder().decode(Wire.self, from: Data(json.utf8)) else {
      return bad("could not read the export request")
    }
    var o = MovieExportOptions()
    switch w.format {
    case "gif": o.format = .gif
    case "png": o.format = .png
    case "mp4", "mov": o.format = .video
    default: return bad("unsupported format '\(w.format)'")
    }
    if o.format == .video {
      switch w.codec {
      case "h264": o.codec = .h264
      case "hevc": o.codec = .hevc
      case "prores": o.codec = .prores
      default: return bad("unsupported video codec '\(w.codec)'")
      }
    }
    let maxDim = MovieExportOptions.maxRenderDimension
    guard (16...maxDim).contains(w.width), (16...maxDim).contains(w.height) else {
      return bad("width and height must be between 16 and \(maxDim) on this device")
    }
    guard [1, 2, 4].contains(w.supersample) else { return bad("supersample must be 1, 2 or 4") }
    guard w.bitrate.isFinite, w.bitrate == 0 || (1...400).contains(w.bitrate) else {
      return bad("bitrate must be 0 (automatic) or 1–400 Mbit/s")
    }
    o.movContainer = (w.format == "mov")
    o.width = w.width
    o.height = w.height
    o.quality = MovieQuality(rawValue: w.quality) ?? .custom
    o.overrides = w.overrides
    // Over the render limit at this size → start() renders at 1×.
    o.supersample = w.supersample
    o.fpsOverride = w.fps
    o.bitrateMbps = w.bitrate
    return .success(
      MovieExportRequest(
        path: w.path, options: o, first: w.first,
        last: w.last, frames: w.frames))
  }
}

// MARK: - Exporter

// NOT @MainActor: the per-frame core render (renderHiResPNG) + AV/GIF encode run
// on `renderQueue` so they don't block the UI (#58 L-59). @Published mutations
// are explicitly hopped back to the main thread; loop control (idx/isExporting)
// is only touched on the main thread (in renderNext / start / finish).
final class MovieExporter: ObservableObject {
  // Serial: frames render + encode one at a time, off the main thread. Also
  // the teardown barrier: anything queued after the in-flight frame runs only
  // once that frame's off-main render has finished.
  private let renderQueue = DispatchQueue(label: "com.mmol.movieexport", qos: .userInitiated)

  @Published var isExporting = false
  @Published var progress: Double = 0
  @Published var finishedURL: URL?
  @Published var errorText: String?
  @Published var etaText: String?
  @Published var estimateText: String?
  @Published var isEstimating = false
  // Encoder finalization in progress: too late to cancel.
  @Published var isFinalizing = false

  // Called on the main thread when an export ends: (result, nil) on success,
  // (nil, message) on failure or cancel. Used by cmd.movie_export, which has
  // no sheet to observe finishedURL.
  var onFinish: ((URL?, String?) -> Void)?

  private weak var engine: PyMOLEngine?
  private var options = MovieExportOptions()
  private var width = 1280, height = 720  // output size
  private var renderW = 1280, renderH = 720  // offscreen render size (× supersample)
  private var first = 1, last = 1, fps = 30
  private var idx = 0
  private var frameDir: URL?
  private var outURL: URL?
  // When frame 1 finished. The ETA averages frames 2+ only, so frame 1's
  // one-off rep rebuild isn't counted once per remaining frame.
  private var firstFrameDone: Date?
  // Set by the first complete()/fail(); later calls are ignored, so a late
  // writer callback can't tear down (or report) an export twice.
  private var ended = false
  // True between apply_overrides and restore; drives every teardown path.
  private var overridesApplied = false
  // True while THIS exporter holds engine.exportRenderActive, so only the
  // owner ever clears it (the sheet and cmd.movie_export each have one).
  private var ownsCore = false

  private var writer: AVAssetWriter?
  private var videoInput: AVAssetWriterInput?
  private var adaptor: AVAssetWriterInputPixelBufferAdaptor?
  private var gifDest: CGImageDestination?

  private var total: Int { max(last - first + 1, 1) }

  private static let restorePython =
    "from pymol import movie_exporting as _me\n_me.restore()"

  func start(
    engine: PyMOLEngine, options: MovieExportOptions,
    first: Int, last: Int, fps: Int
  ) {
    guard !isExporting, !isEstimating else { return }
    guard !engine.exportRenderActive else {
      errorText = "Another movie export is already running."
      onFinish?(nil, errorText)
      return
    }
    self.engine = engine
    self.options = options
    // H.264 requires even dimensions.
    self.width = max(2, options.width - (options.width % 2))
    self.height = max(2, options.height - (options.height % 2))
    let ss = options.allows(supersample: options.supersample) ? options.supersample : 1
    self.renderW = width * ss
    self.renderH = height * ss
    self.first = max(1, min(first, last))
    self.last = max(self.first, last)
    self.fps = max(1, fps)
    self.idx = self.first
    self.progress = 0
    self.errorText = nil
    self.etaText = nil
    self.finishedURL = nil
    self.ended = false
    self.isFinalizing = false
    self.firstFrameDone = nil

    let tmp = FileManager.default.temporaryDirectory
    frameDir = tmp.appendingPathComponent("pymol_frames_\(UUID().uuidString.prefix(6))")
    try? FileManager.default.createDirectory(at: frameDir!, withIntermediateDirectories: true)
    if options.format == .png {
      outURL = tmp.appendingPathComponent("mMol_movie_frames")
    } else {
      outURL = tmp.appendingPathComponent("mMol_movie.\(options.fileExtension)")
    }
    try? FileManager.default.removeItem(at: outURL!)

    guard setupEncoder() else {
      fail(
        "Could not initialize the \(options.format == .video ? options.codec.rawValue : options.format.rawValue) encoder."
      )
      return
    }

    engine.pause()  // stop core-driven advance during capture
    // Export-only overrides, set ONCE before frame 1 so the reps rebuild
    // once (in frame 1's on-main updateScene), not per frame.
    applyOverrides(options.overrides, on: engine)
    // Claim the core for the exporter: the live draw loop + feedback poll now
    // skip (they gate on this), so renderHiResPNG can run off-main exclusively.
    engine.exportRenderActive = true
    ownsCore = true
    isExporting = true
    renderNext()
  }

  // Stop after the in-flight frame. Settings are restored and the partial file
  // discarded once that frame's off-main render has finished.
  func cancel() {
    guard isExporting, !isFinalizing else { return }
    fail("Export cancelled.")
  }

  private func applyOverrides(_ overrides: [String: Int], on engine: PyMOLEngine) {
    guard !overrides.isEmpty,
      let data = try? JSONSerialization.data(withJSONObject: overrides, options: [.sortedKeys])
    else { return }
    engine.runPython(
      "from pymol import movie_exporting as _me\n"
        + "_me.apply_overrides(\(pyLiteral(String(decoding: data, as: UTF8.self))))")
    overridesApplied = true
  }

  private func setupEncoder() -> Bool {
    guard let outURL = outURL else { return false }
    switch options.format {
    case .video:
      let isMOV = outURL.pathExtension == "mov"
      guard let w = try? AVAssetWriter(outputURL: outURL, fileType: isMOV ? .mov : .mp4) else {
        return false
      }
      var settings: [String: Any] = [
        AVVideoCodecKey: {
          switch options.codec {
          case .h264: return AVVideoCodecType.h264
          case .hevc: return AVVideoCodecType.hevc
          case .prores: return AVVideoCodecType.proRes422
          }
        }(),
        AVVideoWidthKey: width, AVVideoHeightKey: height,
      ]
      if options.bitrateMbps > 0 && options.codec != .prores {
        settings[AVVideoCompressionPropertiesKey] =
          [AVVideoAverageBitRateKey: Int(options.bitrateMbps * 1_000_000)]
      }
      let input = AVAssetWriterInput(mediaType: .video, outputSettings: settings)
      input.expectsMediaDataInRealTime = false
      let adaptor = AVAssetWriterInputPixelBufferAdaptor(
        assetWriterInput: input,
        sourcePixelBufferAttributes: [
          kCVPixelBufferPixelFormatTypeKey as String: kCVPixelFormatType_32ARGB,
          kCVPixelBufferWidthKey as String: width,
          kCVPixelBufferHeightKey as String: height,
        ])
      guard w.canAdd(input) else { return false }
      w.add(input)
      guard w.startWriting() else { return false }
      w.startSession(atSourceTime: .zero)
      writer = w
      videoInput = input
      self.adaptor = adaptor
      return true
    case .gif:
      guard
        let dest = CGImageDestinationCreateWithURL(
          outURL as CFURL, UTType.gif.identifier as CFString, total, nil)
      else { return false }
      let props = [kCGImagePropertyGIFDictionary: [kCGImagePropertyGIFLoopCount: 0]]
      CGImageDestinationSetProperties(dest, props as CFDictionary)
      gifDest = dest
      return true
    case .png:
      return
        (try? FileManager.default.createDirectory(
          at: outURL, withIntermediateDirectories: true)) != nil
    }
  }

  // Drives one frame per cycle. Loop control (idx/isExporting/progress) is only
  // touched here on the MAIN thread; the heavy per-frame work — set frame +
  // renderHiResPNG (a blocking GPU render that reshapes the core) + encode —
  // runs on renderQueue OFF the main thread, so the UI stays responsive even
  // for slow frames. The live draw loop + feedback poll are gated
  // by engine.exportRenderActive, so the exporter is the sole core user. (#58 L-59)
  private func renderNext() {
    guard isExporting, let engine = engine, let frameDir = frameDir else { return }
    if idx > last {
      finish()
      return
    }
    let captureIdx = idx
    let w = renderW
    let h = renderH
    let png = frameDir.appendingPathComponent("f\(captureIdx).png")
    // Both of these MUST run on the main thread because they reach the Python
    // C-API, and under _PYMOL_EMBEDDED the main thread owns the interpreter's
    // GIL persistently (PAutoBlock is a no-op, NOT PyGILState_Ensure) — driving
    // Python from the render queue corrupts the Python heap (SIGSEGV in
    // _PyObject_Malloc):
    //   1. cmd.frame(N) — advances to this frame's state.
    //   2. updateScene() — rebuilds this frame's dirty reps now, on-main. The
    //      rebuild goes ObjectMolecule::update -> OrthoBusyFast ->
    //      PLockStatusAttempt (Python), so it must NOT happen inside the
    //      off-main render. Doing it here leaves the off-main SceneRenderMetal
    //      with clean reps and no Python touch.
    // The frame is fully set and rebuilt before the render is dispatched, so
    // there's no overlap (the next frame is only set after this render ends).
    engine.runPython("from pymol import cmd as _c\n_c.frame(\(captureIdx))")
    engine.updateScene()
    renderQueue.async { [weak self, weak engine] in
      guard let self = self, let engine = engine else { return }
      // Off main, exclusive core access (live draw loop + feedback poll are
      // gated by exportRenderActive). The reps were already rebuilt on-main
      // (updateScene above), so this render's SceneUpdate is a clean no-op —
      // pure C++/Metal, no Python — safe off the main thread. It blocks on
      // the GPU and writes the PNG while the UI stays responsive.
      // #601: a video / GIF frame's PNG is read straight back and deleted,
      // so skip its compression (~94% of a 4K frame's time); a PNG
      // sequence keeps normal files. renderHiResPNG has written the file
      // when it returns, so the switch only spans this frame's write.
      let throwaway = self.options.format != .png
      if throwaway { PyMOLBridge_SetFastPNG(1) }
      engine.renderHiResPNG(png.path, width: w, height: h)
      if throwaway { PyMOLBridge_SetFastPNG(0) }
      let error = self.appendFrame(png, frameIndex: captureIdx)  // encode off-main
      try? FileManager.default.removeItem(at: png)
      DispatchQueue.main.async {
        guard self.isExporting else { return }
        if let error = error {
          self.fail("Frame \(captureIdx): \(error)")
          return
        }
        let done = captureIdx - self.first + 1
        self.progress = Double(done) / Double(self.total)
        let now = Date()
        if done == 1 { self.firstFrameDone = now }
        let left = self.total - done
        if left > 0, done >= 2, let t1 = self.firstFrameDone {
          let perFrame = now.timeIntervalSince(t1) / Double(done - 1)
          self.etaText = "About \(Self.formatDuration(perFrame * Double(left))) left"
        } else {
          self.etaText = nil
        }
        self.idx += 1
        self.renderNext()
      }
    }
  }

  // Runs on renderQueue (off main). `frameIndex` is passed in rather than read
  // from `self.idx` (which the main thread mutates) so there's no cross-thread
  // read of the loop counter. Returns why the frame couldn't be written, or nil.
  private func appendFrame(_ png: URL, frameIndex: Int) -> String? {
    guard FileManager.default.fileExists(atPath: png.path) else {
      return "the renderer produced no image"
    }
    switch options.format {
    case .video:
      guard let input = videoInput, let adaptor = adaptor else { return "no encoder" }
      guard let cg = loadCGImage(png) else { return "could not read the rendered image" }
      // Appending while not ready raises; wait (bounded) for the encoder.
      var waited = 0
      while !input.isReadyForMoreMediaData && waited < 30_000 {
        usleep(2000)
        waited += 2
      }
      guard input.isReadyForMoreMediaData else { return "the encoder stopped accepting frames" }
      guard let pb = pixelBuffer(from: cg) else { return "could not convert the image" }
      let t = CMTime(value: Int64(frameIndex - first), timescale: Int32(fps))
      guard adaptor.append(pb, withPresentationTime: t) else {
        return writer?.error?.localizedDescription ?? "the encoder rejected the frame"
      }
    case .gif:
      guard let dest = gifDest else { return "no encoder" }
      guard let cg = loadCGImage(png) else { return "could not read the rendered image" }
      // renderHiResPNG produces a transparent background (molecule alpha=1,
      // bg alpha=0). GIF would collapse that to its background color, so
      // flatten onto opaque black first (matching the MP4 pixel-buffer path).
      let opaque = flattenedOpaque(cg) ?? cg
      let props = [
        kCGImagePropertyGIFDictionary:
          [kCGImagePropertyGIFDelayTime: 1.0 / Double(fps)]
      ]
      CGImageDestinationAddImage(dest, opaque, props as CFDictionary)
    case .png:
      guard let dir = outURL else { return "no output folder" }
      let dst = dir.appendingPathComponent(String(format: "frame_%04d.png", frameIndex - first + 1))
      if renderW == width && renderH == height {
        // No supersampling: keep the renderer's PNG (alpha and all) as-is.
        do { try FileManager.default.moveItem(at: png, to: dst) } catch {
          return error.localizedDescription
        }
      } else {
        guard let cg = loadCGImage(png), let small = downsampled(cg) else {
          return "could not downsample the rendered image"
        }
        guard writePNG(small, to: dst) else { return "could not write \(dst.lastPathComponent)" }
      }
    }
    return nil
  }

  // Composite a (possibly transparent) frame onto opaque black at output size.
  private func flattenedOpaque(_ image: CGImage) -> CGImage? {
    guard
      let ctx = CGContext(
        data: nil, width: width, height: height, bitsPerComponent: 8,
        bytesPerRow: 0, space: CGColorSpaceCreateDeviceRGB(),
        bitmapInfo: CGImageAlphaInfo.noneSkipFirst.rawValue)
    else { return nil }
    if renderW != width { ctx.interpolationQuality = .high }
    ctx.setFillColor(red: 0, green: 0, blue: 0, alpha: 1)
    ctx.fill(CGRect(x: 0, y: 0, width: width, height: height))
    ctx.draw(image, in: CGRect(x: 0, y: 0, width: width, height: height))
    return ctx.makeImage()
  }

  // Scale a supersampled frame down to output size, keeping alpha.
  private func downsampled(_ image: CGImage) -> CGImage? {
    guard
      let ctx = CGContext(
        data: nil, width: width, height: height, bitsPerComponent: 8,
        bytesPerRow: 0, space: CGColorSpaceCreateDeviceRGB(),
        bitmapInfo: CGImageAlphaInfo.premultipliedLast.rawValue)
    else { return nil }
    ctx.interpolationQuality = .high
    ctx.draw(image, in: CGRect(x: 0, y: 0, width: width, height: height))
    return ctx.makeImage()
  }

  @discardableResult
  private func writePNG(_ image: CGImage, to url: URL) -> Bool {
    guard
      let dest = CGImageDestinationCreateWithURL(
        url as CFURL, UTType.png.identifier as CFString, 1, nil)
    else { return false }
    CGImageDestinationAddImage(dest, image, nil)
    return CGImageDestinationFinalize(dest)
  }

  private func finish() {
    switch options.format {
    case .video:
      guard let w = writer else {
        fail("Encoding failed.")
        return
      }
      isFinalizing = true  // Cancel is disabled from here on
      videoInput?.markAsFinished()
      w.finishWriting { [weak self] in
        DispatchQueue.main.async {
          guard let self = self, !self.ended else { return }
          if w.status == .completed {
            self.complete(self.outURL)
          } else {
            self.fail(w.error?.localizedDescription ?? "Encoding failed.")
          }
        }
      }
    case .gif:
      if let dest = gifDest, CGImageDestinationFinalize(dest) {
        complete(outURL)
      } else {
        fail("GIF encoding failed.")
      }
    case .png:
      complete(outURL)
    }
  }

  private func complete(_ url: URL?) {
    guard !ended else { return }
    ended = true
    isFinalizing = false
    isExporting = false
    progress = 1
    etaText = nil
    releaseCore { [weak self] in
      guard let self = self else { return }
      self.finishedURL = url
      self.onFinish?(url, nil)
    }
  }

  private func fail(_ message: String) {
    guard !ended else { return }
    ended = true
    isFinalizing = false
    isExporting = false
    etaText = nil
    errorText = message
    // Cancel the writer, drop the partial output and clear the encoder state
    // on renderQueue, AFTER any in-flight frame has finished using them
    // (appendFrame reads writer/videoInput/adaptor/gifDest off-main).
    let out = outURL
    renderQueue.async { [weak self] in
      if let w = self?.writer, w.status == .writing { w.cancelWriting() }
      if let out = out { try? FileManager.default.removeItem(at: out) }
      self?.writer = nil
      self?.videoInput = nil
      self?.adaptor = nil
      self?.gifDest = nil
    }
    releaseCore { [weak self] in self?.onFinish?(nil, message) }
  }

  // Hand the core back to the live loop: restore the overridden settings and
  // clear exportRenderActive. Queued behind any in-flight frame on renderQueue,
  // because restoring marks reps dirty — if that happened while an off-main
  // render was still running, its SceneUpdate would rebuild reps (and touch
  // Python) off the main thread.
  private func releaseCore(then done: @escaping () -> Void = {}) {
    let engine = self.engine
    let restore = overridesApplied
    let release = ownsCore
    overridesApplied = false
    ownsCore = false
    let dir = frameDir
    frameDir = nil
    renderQueue.async {
      DispatchQueue.main.async {
        if restore { engine?.runPython(MovieExporter.restorePython) }
        if release { engine?.exportRenderActive = false }
        if let d = dir { try? FileManager.default.removeItem(at: d) }
        done()
      }
    }
  }

  // Safety net: if the sheet is dismissed mid-export and the exporter is
  // released, never leave the core flag stuck true (that would freeze the live
  // viewport) or the session's settings overridden. Both wait for the in-flight
  // renderQueue frame, which finishes harmlessly.
  deinit {
    guard let engine = engine, ownsCore || overridesApplied else { return }
    let restore = overridesApplied
    let release = ownsCore
    renderQueue.async {
      DispatchQueue.main.async {
        if restore { engine.runPython(MovieExporter.restorePython) }
        if release { engine.exportRenderActive = false }
      }
    }
  }

  // MARK: time estimate

  // Render one probe frame of the CURRENT frame with these options and push
  // it through the selected encoder/writer (a throwaway one-frame file), then
  // estimate the whole export: the one-off rep rebuild the overrides cause,
  // plus frames × (render + encode/write). Overrides are restored after.
  func estimate(engine: PyMOLEngine, options: MovieExportOptions, frames: Int) {
    guard !isExporting, !isEstimating, !engine.exportRenderActive else { return }
    self.engine = engine
    isEstimating = true
    estimateText = nil
    let w = max(2, options.width - (options.width % 2))
    let h = max(2, options.height - (options.height % 2))
    let ss = options.allows(supersample: options.supersample) ? options.supersample : 1
    self.width = w
    self.height = h
    self.renderW = w * ss
    self.renderH = h * ss
    self.options = options
    self.first = 1
    self.last = 1
    self.fps = max(options.fpsOverride, 30)
    let tmp = FileManager.default.temporaryDirectory
    let tag = UUID().uuidString.prefix(6)
    let png = tmp.appendingPathComponent("pymol_probe_\(tag).png")
    outURL =
      options.format == .png
      ? tmp.appendingPathComponent("pymol_probe_\(tag)_frames")
      : tmp.appendingPathComponent("pymol_probe_\(tag).\(options.fileExtension)")
    let probeOut = outURL
    guard setupEncoder() else {
      isEstimating = false
      estimateText = "Couldn't start the encoder for an estimate."
      return
    }

    engine.pause()
    let t0 = Date()
    applyOverrides(options.overrides, on: engine)
    engine.updateScene()
    let rebuild = Date().timeIntervalSince(t0)
    engine.exportRenderActive = true
    ownsCore = true
    renderQueue.async { [weak self] in
      guard let self = self else { return }
      let t1 = Date()
      engine.renderHiResPNG(png.path, width: w * ss, height: h * ss)
      let error = self.appendFrame(png, frameIndex: 1)
      // Drain the encoder so its compression time is counted too.
      switch options.format {
      case .video:
        let done = DispatchSemaphore(value: 0)
        self.videoInput?.markAsFinished()
        if let wr = self.writer {
          wr.finishWriting { done.signal() }
          done.wait()
        }
      case .gif:
        if let d = self.gifDest { CGImageDestinationFinalize(d) }
      case .png:
        break
      }
      let perFrame = Date().timeIntervalSince(t1)
      try? FileManager.default.removeItem(at: png)
      if let o = probeOut { try? FileManager.default.removeItem(at: o) }
      self.writer = nil
      self.videoInput = nil
      self.adaptor = nil
      self.gifDest = nil
      DispatchQueue.main.async {
        let total = rebuild + perFrame * Double(max(frames, 1))
        self.releaseCore {
          self.isEstimating = false
          self.estimateText =
            error.map { "Estimate failed: \($0)" }
            ?? String(
              format: "≈ %.1f s/frame · about %@ total",
              perFrame, Self.formatDuration(total))
        }
      }
    }
  }

  static func formatDuration(_ s: Double) -> String {
    if s < 60 { return "\(max(1, Int(s.rounded()))) s" }
    let m = Int((s / 60).rounded())
    if m < 60 { return "\(m) min" }
    return String(format: "%d h %02d min", m / 60, m % 60)
  }

  // MARK: helpers

  private func loadCGImage(_ url: URL) -> CGImage? {
    guard let src = CGImageSourceCreateWithURL(url as CFURL, nil) else { return nil }
    return CGImageSourceCreateImageAtIndex(src, 0, nil)
  }

  private func pixelBuffer(from image: CGImage) -> CVPixelBuffer? {
    let attrs: CFDictionary =
      [
        kCVPixelBufferCGImageCompatibilityKey: true,
        kCVPixelBufferCGBitmapContextCompatibilityKey: true,
      ] as CFDictionary
    var pb: CVPixelBuffer?
    CVPixelBufferCreate(kCFAllocatorDefault, width, height, kCVPixelFormatType_32ARGB, attrs, &pb)
    guard let buffer = pb else { return nil }
    CVPixelBufferLockBaseAddress(buffer, [])
    defer { CVPixelBufferUnlockBaseAddress(buffer, []) }
    guard
      let ctx = CGContext(
        data: CVPixelBufferGetBaseAddress(buffer),
        width: width, height: height, bitsPerComponent: 8,
        bytesPerRow: CVPixelBufferGetBytesPerRow(buffer),
        space: CGColorSpaceCreateDeviceRGB(),
        bitmapInfo: CGImageAlphaInfo.noneSkipFirst.rawValue)
    else { return nil }
    // Downsampling a supersampled frame. Left at the default for 1× so a
    // Standard export stays byte-identical to the pre-#581 output.
    if renderW != width { ctx.interpolationQuality = .high }
    ctx.draw(image, in: CGRect(x: 0, y: 0, width: width, height: height))
    return buffer
  }
}

// MARK: - Reusable controls

// The export form, in the Movie tab's right column. Self-contained (owns its
// MovieExporter). Renders the full timeline — no frame-range picker.
struct MovieExportControls: View {
  @EnvironmentObject var engine: PyMOLEngine
  @StateObject private var exporter = MovieExporter()

  private struct SizePreset: Identifiable {
    let name: String
    let w: Int
    let h: Int
    var id: String { name }
  }
  private static let presets: [SizePreset] = [
    SizePreset(name: "360p", w: 640, h: 360),
    SizePreset(name: "480p", w: 854, h: 480),
    SizePreset(name: "720p", w: 1280, h: 720),
    SizePreset(name: "1080p", w: 1920, h: 1080),
    SizePreset(name: "1440p", w: 2560, h: 1440),
    SizePreset(name: "4K", w: 3840, h: 2160),
  ]
  private static let customSizeTag = "custom"

  // Last-used settings, per user (not in the .pse).
  @AppStorage("mmol.movieExport.options") private var storedOptions = Data()
  @State private var options = MovieExportOptions()
  @State private var sizeTag = "720p"
  @State private var customW = 1920
  @State private var customH = 1080
  @State private var showAdvanced = false
  @State private var advancedHeight: CGFloat = 0
  @State private var loaded = false

  private var frameCount: Int { max(engine.playback.frameCount, 1) }
  private var timelineFPS: Int { Int(engine.playback.movieFPS.rounded()) }
  private var effectiveFPS: Int { options.fpsOverride > 0 ? options.fpsOverride : timelineFPS }

  var body: some View {
    VStack(alignment: .leading, spacing: 0) {
      Text("Export Movie").font(.headline).padding(16)

      // The main choices stay put; only the Advanced knobs scroll.
      mainControls.padding(.horizontal, 16)
      advancedDisclosure
        .padding(.horizontal, 16).padding(.top, 18).padding(.bottom, 16)

      Divider()
      footer.padding(16)
    }
    .onAppear(perform: loadStoredOptions)
    .onChange(of: options) { _ in saveOptions() }
    .onChange(of: sizeTag) { _ in applySize() }
    .onChange(of: customW) { _ in applySize() }
    .onChange(of: customH) { _ in applySize() }
    .onChange(of: exporter.finishedURL) { url in
      if let url = url { deliver(url) }
    }
  }

  // Format, Size and the Quality preset.
  private var mainControls: some View {
    VStack(alignment: .leading, spacing: 18) {
      labeled("Format") {
        ClassicSegmented(
          selection: $options.format,
          options: MovieExportOptions.Format.allCases.map { ($0, $0.rawValue) })
      }
      labeled("Size") {
        HStack {
          ClassicPicker(
            selection: $sizeTag,
            options: Self.presets.map {
              ($0.name, "\($0.name)  ·  \(String($0.w))×\(String($0.h))")
            }
              + [(Self.customSizeTag, "Custom…")])
          if sizeTag == Self.customSizeTag {
            dimensionField($customW)
            Text("×").foregroundStyle(.secondary)
            dimensionField($customH)
          }
          Spacer(minLength: 0)
        }
      }
      VStack(alignment: .leading, spacing: 6) {
        labeled("Quality") {
          ClassicSegmented(
            selection: Binding(
              get: { options.quality },
              set: { options.applyPreset($0) }),
            options: MovieQuality.allCases.map { ($0, $0.label) })
        }
        Text(qualityBlurb).font(.caption).foregroundStyle(.secondary)
      }
    }
  }

  // The individual knobs behind the presets. Editing one switches Quality to
  // Custom. On the Mac they scroll in their own area, as tall as the knobs
  // but capped so the sheet fits on screen (#604, #607).
  private var advancedDisclosure: some View {
    DisclosureGroup(isExpanded: $showAdvanced.animation(.easeInOut(duration: 0.2))) {
      ScrollView {
        advancedControls
          .padding(.top, 10).padding(.trailing, 14)
          .background(
            GeometryReader { g in
              Color.clear.preference(key: AdvancedHeightKey.self, value: g.size.height)
            })
      }
      .scrollIndicators(.visible)  // hints there's more below the cap
      .frame(height: min(advancedHeight, Self.advancedMaxHeight))
      .onPreferenceChange(AdvancedHeightKey.self) { advancedHeight = $0 }
    } label: {
      Text("Advanced")
        .font(.system(size: 13, weight: .medium))
    }
  }

  private struct AdvancedHeightKey: PreferenceKey {
    static var defaultValue: CGFloat = 0
    static func reduce(value: inout CGFloat, nextValue: () -> CGFloat) {
      value = max(value, nextValue())
    }
  }

  // Leaves room for the header, main controls and footer (~560 pt) on the
  // current screen.
  private static var advancedMaxHeight: CGFloat {
    let screen = NSScreen.main?.visibleFrame.height ?? 800
    return min(max(screen - 600, 180), 420)
  }

  // Pinned below the scrolling content: summary, warnings, progress and the
  // Render & Export button (#606).
  private var footer: some View {
    VStack(alignment: .leading, spacing: 10) {
      Text("\(frameCount) frames at \(effectiveFPS) fps.")
        .font(.caption).foregroundStyle(.secondary)
      ForEach(warnings, id: \.self) { w in
        Text(w).font(.caption).foregroundStyle(.orange)
          .fixedSize(horizontal: false, vertical: true)
      }

      if exporter.isExporting {
        VStack(alignment: .leading, spacing: 6) {
          ProgressView(value: exporter.progress)
            .tint(TimelineTheme.accent)
          Text(
            "Rendering frame \(Int(exporter.progress * Double(frameCount)))/\(frameCount)…"
              + (exporter.etaText.map { "  \($0)." } ?? "")
          )
          .font(.caption).foregroundStyle(.secondary)
        }
      }
      if let err = exporter.errorText {
        Text(err).font(.caption).foregroundStyle(.red)
          .fixedSize(horizontal: false, vertical: true)
      }

      HStack(spacing: 10) {
        Button(action: runExport) {
          Label(
            exporter.isExporting ? "Rendering…" : "Render & Export",
            systemImage: "film"
          )
          .font(.system(size: 14, weight: .semibold))
          .frame(maxWidth: .infinity)
          .padding(.vertical, 10)
        }
        .buttonStyle(ClassicButtonStyle(prominent: true))
        .disabled(
          exporter.isExporting || exporter.isEstimating
            || engine.playback.frameCount <= 1)
        if exporter.isExporting {
          Button(exporter.isFinalizing ? "Finishing…" : "Cancel", role: .cancel) {
            exporter.cancel()
          }
          .disabled(exporter.isFinalizing)
        }
      }
    }
  }

  // MARK: advanced

  @ViewBuilder
  private var advancedControls: some View {
    VStack(alignment: .leading, spacing: 14) {
      section("Geometry") {
        knob("Cartoon sampling", "cartoon_sampling", [7, 10, 14, 20])
        knob("Sphere quality", "sphere_quality", [1, 2, 3, 4])
        knob("Stick quality", "stick_quality", [8, 12, 16, 24])
        knob("Surface quality", "surface_quality", [0, 1, 2, 3])
      }
      section("Image") {
        row("Supersampling") {
          ClassicSegmented(
            selection: Binding(
              get: { options.supersample },
              set: {
                options.supersample = $0
                options.quality = .custom
              }),
            options: [1, 2, 4].map { ($0, "\($0)×") }
          )
          .frame(width: 140)
        }
        knob("Anti-aliasing (MSAA)", "metal_msaa", [0, 1], onOff: true)
      }
      if options.format == .video {
        section("Encoding") {
          row("Codec") {
            ClassicPicker(
              selection: $options.codec,
              options: MovieExportOptions.Codec.allCases.map { ($0, $0.rawValue) })
          }
          if options.codec != .prores {
            row("Bitrate") {
              ClassicPicker(
                selection: $options.bitrateMbps,
                options: [(0.0, "Automatic")]
                  + [8.0, 16, 32, 64, 128].map { ($0, "\(Int($0)) Mbit/s") })
            }
          }
        }
      }
      if options.format != .png {
        row("Frame rate") {
          ClassicPicker(
            selection: $options.fpsOverride,
            options: [(0, "Timeline (\(timelineFPS) fps)")]
              + [24, 25, 30, 50, 60].map { ($0, "\($0) fps") })
        }
      }

      HStack(spacing: 8) {
        Button(exporter.isEstimating ? "Estimating…" : "Estimate time") {
          exporter.estimate(engine: engine, options: options, frames: frameCount)
        }
        .disabled(exporter.isExporting || exporter.isEstimating)
        if let e = exporter.estimateText {
          Text(e).font(.caption).foregroundStyle(.secondary)
        }
      }
    }
  }

  private var qualityBlurb: String {
    switch options.quality {
    case .draft: return "Faster: no MSAA. For checking timing and framing."
    case .standard: return "The same settings as the live view."
    case .high: return "Finer geometry."
    case .maximum: return "Finest geometry and surfaces, 2× supersampling. Slow."
    case .custom: return "Your own mix. \"Session\" keeps the current setting."
    }
  }

  private var warnings: [String] {
    var w: [String] = []
    if (options.overrides["surface_quality"] ?? 0) >= 2 {
      w.append(
        "High surface quality rebuilds surfaces at export start. This can take minutes and a lot of memory."
      )
    }
    if !options.allows(supersample: options.supersample) {
      w.append(
        "\(options.supersample)× supersampling exceeds the \(MovieExportOptions.maxRenderDimension)-pixel render limit at this size; frames render at 1×."
      )
    }
    return w
  }

  // A menu of fixed values plus "Session" (nil = don't override). Picking a
  // value switches the preset to Custom.
  private func knob(
    _ title: String, _ key: String, _ values: [Int],
    onOff: Bool = false
  ) -> some View {
    row(title) {
      ClassicPicker(
        selection: Binding<Int?>(
          get: { options.overrides[key] },
          set: { v in
            options.overrides[key] = v
            options.quality = .custom
          }),
        options: [(Int?.none, "Session")]
          + values.map { (Int?.some($0), onOff ? ($0 != 0 ? "On" : "Off") : "\($0)") })
    }
  }

  private func row<C: View>(_ title: String, @ViewBuilder _ content: () -> C) -> some View {
    HStack {
      Text(title).font(.system(size: 12))
      Spacer()
      content()
    }
  }

  private func section<C: View>(_ title: String, @ViewBuilder _ content: () -> C) -> some View {
    VStack(alignment: .leading, spacing: 6) {
      Text(title.uppercased())
        .font(.system(size: 10, weight: .semibold)).foregroundStyle(.secondary)
      content()
    }
  }

  private func dimensionField(_ value: Binding<Int>) -> some View {
    TextField("", value: value, format: .number.grouping(.never))
      .classicField()
      .frame(width: 64)
  }

  // MARK: state

  private func applySize() {
    guard loaded else { return }
    if sizeTag == Self.customSizeTag {
      options.width = min(max(customW, 16), MovieExportOptions.maxRenderDimension)
      options.height = min(max(customH, 16), MovieExportOptions.maxRenderDimension)
    } else if let p = Self.presets.first(where: { $0.name == sizeTag }) {
      options.width = p.w
      options.height = p.h
    }
  }

  private func loadStoredOptions() {
    if let o = try? JSONDecoder().decode(MovieExportOptions.self, from: storedOptions) {
      options = o
    }
    if let p = Self.presets.first(where: { $0.w == options.width && $0.h == options.height }) {
      sizeTag = p.name
    } else {
      sizeTag = Self.customSizeTag
      customW = options.width
      customH = options.height
    }
    // Applied after the size tag so the restored size isn't overwritten.
    DispatchQueue.main.async { loaded = true }
  }

  private func saveOptions() {
    guard loaded, let d = try? JSONEncoder().encode(options) else { return }
    storedOptions = d
  }

  private func runExport() {
    exporter.start(
      engine: engine, options: options,
      first: 1, last: frameCount, fps: effectiveFPS)
  }

  // MARK: deliver result

  private func deliver(_ url: URL) {
    let panel = NSSavePanel()
    panel.nameFieldStringValue = url.lastPathComponent
    if !url.pathExtension.isEmpty, let ct = UTType(filenameExtension: url.pathExtension) {
      panel.allowedContentTypes = [ct]
    }
    panel.canCreateDirectories = true
    panel.title = options.format == .png ? "Save Frames Folder" : "Save Movie"
    guard panel.runModal() == .OK, let dest = panel.url else { return }
    try? FileManager.default.removeItem(at: dest)
    try? FileManager.default.copyItem(at: url, to: dest)
  }

  @ViewBuilder
  private func labeled<C: View>(_ title: String, @ViewBuilder _ content: () -> C) -> some View {
    VStack(alignment: .leading, spacing: 6) {
      Text(title).font(.system(size: 12, weight: .medium)).foregroundStyle(.secondary)
      content()
    }
  }
}
