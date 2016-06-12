#ifndef DOMAIN_H
#define DOMAIN_H

namespace myOctree {

enum domain_flags {

	CORNER,		//used for corners
	EAST_GHOST,
	WEST_GHOST,
	NORTH_GHOST,
	SOUTH_GHOST,
	TOP_GHOST,
	BOTTOM_GHOST,	
	DOMAIN

};

typedef domain_flags Mask;


}
#endif
