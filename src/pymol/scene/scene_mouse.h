#pragma once

#include "pymol/scene/picking.h"
#include "pymol/zstring_view.h"

struct Block;
struct PyMOLGlobals;

namespace pymol {
struct CObject;
}

/**
 * Handles mouse clicks in the scene
 * @param block Scene block
 * @param button button pressed
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneClick(Block* block, int button, int x, int y, int mod, double when);

/**
 * Handles mouse release in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneRelease(Block* block, int button, int x, int y, int mod, double when);

/**
 * Handles mouse drag in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneDrag(Block* block, int x, int y, int mod, double when);

/**
 * Handles mouse drag in the scene
 * @param block Scene block
 * @param button button released
 * @param x x coordinate
 * @param y y coordinate
 * @param mod modifier key
 * @param when time of event
 *
 * @pre Active rendering context
 */
void SceneMouseMove(Block* block, int x, int y, int mod, double when);

/**
 * Deselect or toggle off entries
 * @param button button pressed
 * @param mod modifier key
 * @param mode button mode
 * @return if something has been deselected or toggled
 */

void SceneClickPickNothing(PyMOLGlobals* G, int button, int mod, int mode);
