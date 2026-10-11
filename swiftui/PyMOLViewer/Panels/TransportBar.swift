// TransportBar.swift — the unified Timeline transport for states, NMR models,
// trajectory frames and movie frames (all one cmd.frame index in PyMOL).
//
// Familiar desktop VCR transport (|< < play/pause > >|), a frame counter, a
// loop toggle and the frame rate, embedded in the TimelinePanel (whose ruler is
// the scrubber). Playback is core-driven (engine.play()/pause() →
// cmd.mplay/mstop, advanced by the renderer idle tick); the bar only reflects
// state.

import SwiftUI

enum TimelineTheme {
  // Transport / movie / scene controls in the classic chrome: the accent,
  // the panel grey (slightly translucent, since the docked timeline floats
  // over the viewport) and the panel text at full and reduced strength.
  static let accent = ClassicUI.accent
  static let bar = ClassicUI.panelBg.opacity(0.96)
  static let text = ClassicUI.chromeText.color
  static let dim = ClassicUI.chromeText.color.opacity(0.6)
  static let subtleFill = ClassicUI.chromeText.color.opacity(0.12)
}

// Transport control sizing.
private let kTBtnW: CGFloat = 30
private let kTPlayW: CGFloat = 34
private let kTBtnH: CGFloat = 30

struct TransportBar: View {
  @EnvironmentObject var engine: PyMOLEngine
  // Observe the isolated playback object so frame ticks re-render ONLY this
  // bar (not the inspector/menus, which observe `engine`).
  @EnvironmentObject var playback: PlaybackState

  /// When set (movie/timeline transport), the counter + scrubber reflect the
  /// AUTHORED movie length instead of the core frame count — so the movie
  /// transport stays decoupled from multi-state model inspection.
  var movieFrames: Int? = nil

  var body: some View {
    regularRow.background(TimelineTheme.bar)
  }

  // MARK: - Layouts

  // The ruler is the scrubber, so transport, counter, loop and frame rate
  // share one line.
  private var regularRow: some View {
    HStack(spacing: 10) {
      transportCluster
      counter
      loopButton
      fpsMenu
    }
    .padding(.horizontal, 12)
    .frame(height: 40)
  }

  // MARK: - Pieces

  private var transportCluster: some View {
    HStack(spacing: 2) {
      iconButton("backward.end.fill", size: 15) { engine.rewindMovie() }
        .accessibilityLabel("Rewind to start").help("Rewind to start")
      iconButton("backward.fill", size: 15) { engine.stepBackward() }
        .accessibilityLabel("Step back").help("Step back one frame")
      playPauseButton(size: 20)
      iconButton("forward.fill", size: 15) { engine.stepForward() }
        .accessibilityLabel("Step forward").help("Step forward one frame")
      iconButton("forward.end.fill", size: 15) { engine.endingMovie() }
        .accessibilityLabel("Go to end").help("Jump to end")
    }
  }

  private func playPauseButton(size: CGFloat) -> some View {
    Button(action: { engine.togglePlay() }) {
      Image(systemName: playback.isPlaying ? "pause.fill" : "play.fill")
        .font(.system(size: size, weight: .semibold))
        .foregroundColor(TimelineTheme.accent)
        .frame(width: kTPlayW, height: kTBtnH)
        .contentShape(Rectangle())
    }
    .buttonStyle(.plain)
    .accessibilityLabel(playback.isPlaying ? "Pause" : "Play")
    .help(playback.isPlaying ? "Pause" : "Play")
  }

  // Total frames shown/scrubbed: the authored movie length when provided
  // (timeline transport), else the core frame count.
  private var totalFrames: Int { max(movieFrames ?? playback.frameCount, 1) }
  private var shownFrame: Int { min(max(playback.currentFrame, 1), totalFrames) }

  private var counter: some View {
    // Reserve width for the movie's max digit count ("NNN / NNN") so the
    // current frame number growing (1 → N digits) never shifts neighbors.
    let digits = String(totalFrames).count
    return Text("\(shownFrame) / \(totalFrames)")
      .font(.system(size: 12, weight: .medium, design: .monospaced))
      .foregroundColor(TimelineTheme.text)
      .lineLimit(1)
      .frame(minWidth: CGFloat(digits * 2 + 3) * 7.4, alignment: .trailing)
  }

  private var loopButton: some View {
    Button(action: { engine.setMovieLoop(!playback.movieLoop) }) {
      Image(systemName: "repeat")
        .font(.system(size: 15, weight: .semibold))
        .foregroundColor(playback.movieLoop ? TimelineTheme.accent : TimelineTheme.dim)
        .frame(width: kTBtnW, height: kTBtnH)
        .contentShape(Rectangle())
    }
    .buttonStyle(.plain)
    .accessibilityLabel(playback.movieLoop ? "Looping on" : "Looping off")
    .help(playback.movieLoop ? "Looping on — tap to turn off" : "Loop the movie")
  }

  // The frame-rate choices + the "Show frame rate" HUD toggle (issue #142).
  private var fpsMenu: some View {
    ClassicMenuTrigger {
      Label("\(fpsLabel) fps", systemImage: "gauge.with.dots.needle.67percent")
        .font(.system(size: 12))
        .foregroundColor(TimelineTheme.accent)
    } items: {
      [30.0, 15.0, 5.0, 1.0, 0.3].map { f in
        .check(f == f.rounded() ? "\(Int(f)) fps" : "\(f) fps", playback.movieFPS == f) {
          engine.setMovieFPS(f)
        }
      } + [.separator, .action("Show frame rate") { engine.setShowFrameRate(true) }]
    }
    .help("Playback frame rate")
  }

  private var fpsLabel: String {
    let f = playback.movieFPS
    return f == f.rounded() ? String(Int(f)) : String(format: "%.1f", f)
  }

  private func iconButton(
    _ systemName: String, size: CGFloat,
    action: @escaping () -> Void
  ) -> some View {
    Button(action: action) {
      Image(systemName: systemName)
        .font(.system(size: size, weight: .semibold))
        .foregroundColor(TimelineTheme.text)
        .frame(width: kTBtnW, height: kTBtnH)
        .contentShape(Rectangle())
    }
    .buttonStyle(.plain)
  }
}
