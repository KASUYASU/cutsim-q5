/*
 *  Copyright 2015      Kazuyasu Hamada (k-hamada@gifu-u.ac.jp)
*/

#ifndef DEFINITION_H
#define DEFINITION_H

#define MULTI_AXIS
#define PI 		(3.1415926535897932)
#define CALC_TOLERANCE	(1e-6)

#define X_AXIS		(1)
#define Y_AXIS		(2)
#define Z_AXIS		(3)
#ifdef MULTI_AXIS
#define A_AXIS		(4)
#define B_AXIS		(5)
#define C_AXIS		(6)
#endif

//#define SIGN	(-1.0)
#define SIGN_A	(-1.0)
#define SIGN_B	(-1.0)
#define SIGN_C	(-1.0)

#define TOLERANCE			(1e-2)
#define COLLISION_TOLERANCE	(2e-2)

// Node pooling
#define POOL_NODE

// Multi Threadings
#define MULTI_THREAD_SUM
#define MULTI_THREAD_INTERSECT
//#define MULTI_THREAD_DIFF
#define MULTI_THREAD_STL_NEIGHBOR

#if (defined(MULTI_THREAD_SUM) || defined(MULTI_THREAD_INTERSECT) || defined(MULTI_THREAD_DIFF) || defined(MULTI_THREAD_STL_NEIGHBOR))
	#define MULTI_THREAD
#endif

// Display by Wire frame
//#define WIRE_FRAME

#define DEFAULT_SCENE_RADIUS	(100)

#define DEFAULT_CUBE_SIZE		(100.0)
#define DEFAULT_MAX_DEPTH		(9)

#define DEFAULT_STEP_SIZE		(0.1)

#define TOOL_BODY_COLOR		0.9, 0.9, 0.85
#define TOOL_FLUTE_COLOR		0.7, 0.7, 0.65
#define STOCK_COLOR			0,1,1
// PARTS must be gray color.
#define PARTS_COLOR			1,1,1
#define CUTTING_COLOR			1,1,0
#define COLLISION_COLOR		1,0,0

#define DEFAULT_TRAVERSE_FEED_RATE	(1000.0)
#define DEFAULT_FEED_RATE          	(200.0)

#define DEFAULT_HOLDER_RADIUS		(15.0)
#define DEFAULT_HOLDER_LENGTH		(20.0)
#define DEFAULT_HOLDER_COLOR		0.2, 0.2, 0.2

#define DEFAULT_SPINDLE_RADIUS		(11.0)
#define DEFAULT_SPINDLE_LENGTH		(35.0)
#define DEFAULT_SPINDLE_COLOR		0.9, 0.9, 0.9

#define DEFAULT_MAX_X_LIMIT			(150)
#define DEFAULT_MIN_X_LIMIT			(-150)
#define DEFAULT_MAX_Y_LIMIT			(200)
#define DEFAULT_MIN_Y_LIMIT			(-200)
#define DEFAULT_MAX_Z_LIMIT			(150)
#define DEFAULT_MIN_Z_LIMIT			(0)
#define DEFAULT_MAX_FEED_RATE 		(1500.0)
#define DEFAULT_MAX_SPINDLE_POWER 	(150.0)

#define DEFAULT_MAX_TOOL_SLOT		(100)

#define DEFAULT_SPECIFIC_CUTTING_FORCE	(2000)

#define DEFAULT_ANIMATE_INTERVAL	(3)

#endif // DEFINITION_H
