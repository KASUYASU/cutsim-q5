/*
 *  Copyright 2010-2011 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 *  Copyright 2015      Kazuyasu Hamada (k-hamada@gifu-u.ac.jp)
 *
 *  This file is part of OpenCAMlib.
 *
 *  OpenCAMlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  OpenCAMlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with OpenCAMlib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <list>
#include <cassert>
#include <iostream>
#include <sstream>

#include <boost/foreach.hpp>

#include <src/app/cutsim_def.hpp>

#include "octnode.hpp"
#include "octree.hpp"
#include "volume.hpp"

namespace cutsim {

//**************** Octree ********************/

Octree::Octree(double scale, unsigned int  depth, GLVertex* centerp, GLData* gl) {
    root_scale = scale;
    max_depth = depth;
    g = gl;
    // parent(=root) scale, GLdata, center
    root = new Octnode( centerp, root_scale, g );

    for ( int n=0;n<8;++n) {
        root->child[n] = NULL;
    }
    debug = false;
    debug_mc = false;
for (unsigned int i = 0; i < max_depth; i++) {
	int amount = 0;
	for (int k = (max_depth - 1) - i; k >= 0; k--)
		amount += (int)pow((1 << k), 3);
	processed.push_back(amount);
}
}

Octree::~Octree() {
    delete root;
    root = 0;
}

unsigned int Octree::get_max_depth() const {
    return max_depth;
}

double Octree::get_root_scale() const {
    return root_scale;
}

double Octree::leaf_scale() const {
    return (2.0*root_scale) / pow(2.0, (int)max_depth);
}

/// subdivide the Octree n times
void Octree::init(const unsigned int n) {
    for (unsigned int m=0;m<n;++m) {
        std::vector<Octnode*> nodelist;
        get_leaf_nodes(root, nodelist);
        BOOST_FOREACH( Octnode* node, nodelist) {
            node->force_subdivide();
        }
    }
}

void Octree::get_invalid_leaf_nodes( std::vector<Octnode*>& nodelist) const {
    get_invalid_leaf_nodes( root, nodelist );
}

void Octree::get_invalid_leaf_nodes(Octnode* current, std::vector<Octnode*>& nodelist) const {
    if ( current->childcount == 0 ) {
        if ( !current->valid() ) {
            nodelist.push_back( current );
        }
    } else {//surface()surface()
        for ( int n=0;n<8;++n) {
            if ( current->hasChild(n) ) {
                if ( !current->valid() ) {
                    get_leaf_nodes( current->child[n], nodelist );
                }
            }
        }
    }
}

/// put leaf nodes into nodelist
void Octree::get_leaf_nodes(Octnode* current, std::vector<Octnode*>& nodelist) const {
    if ( current->isLeaf() ) {
        nodelist.push_back( current );
    } else {
        for ( int n=0;n<8;++n) {
            if ( current->child[n] != 0 )
                get_leaf_nodes( current->child[n], nodelist );
        }
    }
}

/// put all nodes into nodelist
void Octree::get_all_nodes(Octnode* current, std::vector<Octnode*>& nodelist) const {
    if ( current ) {
        nodelist.push_back( current );
        for ( int n=0;n<8;++n) {
            if ( current->child[n] != 0 )
                get_all_nodes( current->child[n], nodelist );
        }
    }
}

// sum (union) of tree and Volume
//void Octree::sum(Octnode* current, const Volume* vol) {
void Octree::sum(Octnode* current, Volume* vol) {
//    if ( current->is_inside() || !vol->bb.overlaps( current->bb ) ) // if no overlap, or already INSIDE, then quit.
//        return; // abort if no overlap.
if ( current->is_inside() || !vol->bb.overlaps( current->bb ) ) { // if no overlap, or already INSIDE, then quit.
vol->accumlateProgress(processed[current->depth]);
return; // abort if no overlap.
}

//if ((current->depth == (this->max_depth-1)) && current->is_undecided()) return;
if ((current->depth == (this->max_depth-1)) && current->is_undecided() && !current->color.compareColor(vol->color)) { // return;
	vol->accumlateProgress(1);
	return;
}

    current->sum(vol);

if (vol->type == cutsim::STL_VOLUME) {
    if (current->depth == (this->max_depth-1)) {
		current->set_state();
		vol->accumlateProgress(1);
		return;
	} else if (current->check_complete_inside_outside()) { //return;
		vol->accumlateProgress(processed[current->depth]);
		return;
	}
} else
	 current->set_state();

    if ( (current->childcount == 8) ) { // recurse into existing tree
#ifdef MULTI_THREAD_SUM
            QFuture<void> future[8];
            bool dispatch[8];
            for(int m=0;m<8;++m)
            	if ( !current->child[m]->is_inside()  ) { // nodes that are already INSIDE cannot change in a sum-operation
            		future[m] = QtConcurrent::run(this, &Octree::sum, current->child[m], vol);
            		dispatch[m] = true;
            	} else
            		dispatch[m] = false;
            for(int m=0;m<8;++m)
            	if (dispatch[m] == true)
            		future[m].waitForFinished();
#else
        for(int m=0;m<8;++m) {
            if ( !current->child[m]->is_inside()  ) // nodes that are already INSIDE cannot change in a sum-operation
                sum( current->child[m], vol); // call sum on children
        }
#endif
    } else { // no children, subdivide it
        if ( (current->depth < (this->max_depth-1)) ) {
            if (!current->is_undecided()) { current->force_setUndecided(); }
            current->subdivide(); // smash into 8 sub-pieces
#ifdef MULTI_THREAD_SUM
            QFuture<void> future[8];
            for(int m=0;m<8;++m)
            	future[m] = QtConcurrent::run(this, &Octree::sum, current->child[m], vol);
            for(int m=0;m<8;++m)
            	future[m].waitForFinished();
#else
            for(int m=0;m<8;++m)
                sum( current->child[m], vol); // call sum on children
#endif
        }
    }
    // now all children of current have their status set, and we can prune.
    if ( (current->childcount == 8) && ( current->all_child_state(Octnode::INSIDE) || current->all_child_state(Octnode::OUTSIDE) ) ) {
    	if (current->all_child_state(Octnode::INSIDE))
    		current->state = Octnode::INSIDE;
    	else
    		current->state = Octnode::OUTSIDE;
        current->delete_children();
    }

vol->accumlateProgress(1 << ((this->max_depth-1) - current->depth));
if (current->depth < (this->max_depth-3))
vol->sendProgress();

    return;
}

// diff (intersection with volume's compliment) of tree and Volume
void Octree::diff(Octnode* current, const Volume* vol) {
    if ( current->is_outside() || !vol->bb.overlaps( current->bb ) ) // if no overlap, or already OUTSIDE, then return.
    	return;

    current->diff(vol);

	current->set_state();

    if ( ((current->childcount) == 8) /*&& current->is_undecided()*/ ) { // recurse into existing tree
         for(int m=0;m<8;++m) {
            if ( !current->child[m]->is_outside()  ) // nodes that are OUTSIDE don't change
                diff( current->child[m], vol); // call diff on children
        }
    } else { // no children, subdivide it
        if ( (current->depth < (this->max_depth-1)) ) {
            if (!current->is_undecided()) { current->force_setUndecided(); }
            current->subdivide(); // smash into 8 sub-pieces
            for(int m=0;m<8;++m) {
                diff( current->child[m], vol); // call diff on children
            }
        }
    }
    // now all children have their status set, prune.
    if ( (current->childcount == 8) && ( /*current->all_child_state(Octnode::INSIDE) ||*/ current->all_child_state(Octnode::OUTSIDE) ) ) {
    	current->state = Octnode::OUTSIDE;
        current->delete_children();
    }
}

// intersect (intersection) of tree and Volume
void Octree::intersect(Octnode* current, const Volume* vol) {
    if ( current->is_outside() ) // if already OUTSIDE, then return.
        return;   
    
    current->intersect(vol);

	current->set_state();

    if ( ((current->childcount) == 8) /*&& current->is_undecided()*/ ) { // recurse into existing tree
        for(int m=0;m<8;++m) {
            //if ( !current->child[m]->is_outside()  ) // nodes that are OUTSIDE don't change
                intersect( current->child[m], vol); // call intersect on children
        }
//    } else if (  current->is_undecided() ) { // no children, subdivide if undecided 
//    	if (current->childcount != 0) { std::cout << " current->childcount != 0 now:" << current->childcount; return; }
    } else { // no children, subdivide it
    	if (!current->is_undecided()) { current->force_setUndecided(); }
        if ( (current->depth < (this->max_depth-1)) ) {
            current->subdivide(); // smash into 8 sub-pieces
            for(int m=0;m<8;++m) {
                intersect( current->child[m], vol); // call intersect on children
            }
        }
    }
    // now all children have their status set, prune.
    if ( (current->childcount == 8) && ( current->all_child_state(Octnode::INSIDE) || current->all_child_state(Octnode::OUTSIDE) ) ) {
    	if (current->all_child_state(Octnode::INSIDE))
    		current->state = Octnode::INSIDE;
    	else
    		current->state = Octnode::OUTSIDE;
        current->delete_children();
    }
}

// diff (intersection with volume's compliment) of tree and Volume for cuttings
CuttingStatus Octree::diff_c(Octnode* current, const Volume* vol) {
    CuttingStatus status = { 0, NO_COLLISION }, childstatus;
    if ( current->is_outside() || (!vol->bb.overlaps( current->bb ) && (!((CutterVolume*)vol)->enableholder || !((CutterVolume*)vol)->bbHolder.overlaps( current->bb ))) )
    	return status;

    if (current->depth == (this->max_depth-1)) {
    	status = current->diff_cd(vol);
    	if (current->set_state() == false) {
    		status.cutcount = 0;
    		status.collision = NO_COLLISION;
    	}
    	return status;
    } else {
    	current->diff(vol);
    	current->set_state();
    }

    if ( ((current->childcount) == 8) /*&& current->is_undecided()*/ ) { // recurse into existing tree
#ifdef MULTI_THREAD_DIFF
    	QFuture<CuttingStatus> future[8];
        bool dispatch[8];
        for(int m=0;m<8;++m) {
            if ( !current->child[m]->is_outside() ) {
            	future[m] = QtConcurrent::run(this, &Octree::diff_c, current->child[m], vol);
            	dispatch[m] = true;
            } else
            	dispatch[m] = false;
        }
        for(int m=0;m<8;++m) {
        	if (dispatch[m] == true) {
            	future[m].waitForFinished();
            	childstatus = future[m].result();
            	status.cutcount += childstatus.cutcount;
            	status.collision |= childstatus.collision;
            }
        }
#else
        for(int m=0;m<8;++m) {
        	if ( !current->child[m]->is_outside() ) { // nodes that are OUTSIDE don't change
        		childstatus = diff_c( current->child[m], vol); // call diff_c on children
        		status.cutcount += childstatus.cutcount;
        		status.collision |= childstatus.collision;
        	}
        }
#endif
    } else { // no children, subdivide it
         if ( (current->depth < (this->max_depth-1)) ) {
    		if (!current->is_undecided()) { current->force_setUndecided(); }
    		current->subdivide(); // smash into 8 sub-pieces
#ifdef MULTI_THREAD_DIFF
            QFuture<CuttingStatus> future[8];
            for(int m=0;m<8;++m)
            	future[m] = QtConcurrent::run(this, &Octree::diff_c, current->child[m], vol);
            for(int m=0;m<8;++m) {
            	future[m].waitForFinished();
            	childstatus = future[m].result();
				status.cutcount += childstatus.cutcount;
				status.collision |= childstatus.collision;
            }
#else
    		for(int m=0;m<8;++m) {
    			childstatus = diff_c( current->child[m], vol); // call diff_c on children
    			status.cutcount += childstatus.cutcount;
    			status.collision |= childstatus.collision;
    		}
#endif
         }
    }
    // now all children have their status set, prune.
    if ( (current->childcount == 8) && ( /*current->all_child_state(Octnode::INSIDE) ||*/ current->all_child_state(Octnode::OUTSIDE) ) ) {
    	current->state = Octnode::OUTSIDE;
        current->delete_children();
    }
    return status;
}

bool Octree::check_node(Octnode* current, const Volume* vol)
{
	if (vol->type != cutsim::STL_VOLUME)
		return true;

    if ((current->depth == (this->max_depth-1)) && current->is_undecided()) {
    	return current->check_f_value_with_limit();
    }

	if ( ((current->childcount) == 8) && current->is_undecided() ) { // recurse into existing tree
		int state = current->check_node_state();
		if (state == Octnode::UNDECIDED) state = current->parent->check_node_state();
		if (state == Octnode::UNDECIDED) {
			bool outside = true;
			bool inside  = true;
			for(int m=0;m<8;++m) {
				if (current->child[m]->is_inside()) outside = false;
			    if (current->child[m]->is_outside()) inside = false;
			}
			state = inside ? Octnode::INSIDE : outside ? Octnode::OUTSIDE : Octnode::UNDECIDED;
		}

        for(int m=0;m<8;++m) {
        	bool check = check_node(current->child[m], vol); // call check_node on children
        	if (check == false) {
           		if (state == Octnode::OUTSIDE) {
        			std::cout << "Force Outside x:" << current->child[m]->center->x << " y:" << current->child[m]->center->y << " z:" << current->child[m]->center->z << "\n";
        			current->child[m]->color.set(1,1,1);
//        			current->child[m]->state = Octnode::OUTSIDE;
        		} else if (state == Octnode::INSIDE) {
        			std::cout << "Force Inside  x:" << current->child[m]->center->x << " y:" << current->child[m]->center->y << " z:" << current->child[m]->center->z << "\n";
        			current->child[m]->color.set(1,0,0);
//        			current->child[m]->state = Octnode::INSIDE;
        		} else
        			std::cout << "Can't Decide  x:" << current->child[m]->center->x << " y:" << current->child[m]->center->y << " z:" << current->child[m]->center->z << "\n";
        	}
    	}
    }

    // now all children have their status set, prune.
//    if ( (current->childcount == 8) && ( current->all_child_state(Octnode::INSIDE) || current->all_child_state(Octnode::OUTSIDE) ) ) {
    if ( (current->childcount == 8) && (current->check_include_undecided() == false) ) {
    	if (current->all_child_state(Octnode::INSIDE))
    		current->state = Octnode::INSIDE;
//    	else
    	else if (current->all_child_state(Octnode::OUTSIDE))
    		current->state = Octnode::OUTSIDE;
    	else {
    		int inside = 0, outside = 0;
			for(int m=0;m<8;++m) {
				if (current->child[m]->is_inside())  outside++;
			    if (current->child[m]->is_outside()) inside--;
			}
			std::cout << "INSIDE mix OUTSIDE  x:" << current->center->x << " y:" << current->center->y << " z:" << current->center->z << " inside " << inside << " outside " << outside << "\n";
			goto skip_delete_children;
    	}
		current->delete_children();
    }
skip_delete_children:
    return true;
}

#ifdef POOL_NODE
extern std::vector<Octnode*> nodePool;
#endif

// string repr
std::string Octree::str() const {
    std::ostringstream o;
    o << " Octree: ";
    std::vector<Octnode*> nodelist;
    Octree::get_all_nodes(root, nodelist);
    std::vector<int> nodelevel(this->max_depth);
    std::vector<int> invalidsAtLevel(this->max_depth);
    std::vector<int> surfaceAtLevel(this->max_depth);
    int totalVertexSize = 0;
    BOOST_FOREACH( Octnode* n, nodelist) {
        ++nodelevel[n->depth];
        if ( !n->valid() )
            ++invalidsAtLevel[n->depth];
        if (n->is_undecided() )
            ++surfaceAtLevel[n->depth];
        for (int i = 0; i < 8; i++)
        	 if (n->vertex[i] != 0)
        		 totalVertexSize += sizeof(GLVertex);
    }
    o << "  " << nodelist.size() << " leaf-nodes:\n";
    int m=0;
    BOOST_FOREACH( int count, nodelevel) {
        o << "depth="<<m <<"  " << count << " nodes, " << invalidsAtLevel[m] << " invalid, surface=" << surfaceAtLevel[m] << " \n";
        ++m;
    }
    o << "  total " << nodelist.size() * sizeof (Octnode) << " bytes comsumed for node(" << sizeof(Octnode) << " bytes). \n";
    extern unsigned int alocation_count;
    extern unsigned int delete_count;
    extern unsigned int delete_childlen_count;
    o << "    alocation count: " << alocation_count << "  delete count: " << delete_count << "  difference: " << alocation_count - delete_count << "\n";
    o << "    delete child count: " << delete_childlen_count << "\n";
    o << "  total " << totalVertexSize << " bytes comsumed for vertex(" << sizeof(GLVertex) << " bytes). \n";
#ifdef POOL_NODE
    o << "  Node Pool size " << nodePool.size() << "\n";
#endif
    return o.str();
}

void Octree::treeTransfer(Octnode* current, GLVertex parallel, int flip_axis, bool ignore_parts) {
	current->nodeTransfer(parallel, flip_axis, ignore_parts);
	if (((current->childcount) == 8))
		for(int m=0;m<8;++m)
			treeTransfer(current->child[m], parallel, flip_axis, ignore_parts);
}

void Octree::setInvalid(Octnode* current) {
	if (current->is_undecided() && current->isLeaf() && current->valid())
		current->setInvalid();
	if (((current->childcount) == 8))
		for(int m=0;m<8;++m)
			setInvalid(current->child[m]);
}


void Octree::clearTree(double root_scale, unsigned int max_depth, GLVertex* centerPoint) {
	clearNode(root);
    root->scale = root_scale;
    root->depth = 0;

    root->center = centerPoint;
    root->state = Octnode::UNDECIDED;
    root->prev_state = Octnode::OUTSIDE;

    const GLVertex direction[8] = {
                         GLVertex( 1, 1,-1),   // 0
                         GLVertex(-1, 1,-1),   // 1
                         GLVertex(-1,-1,-1),   // 2
                         GLVertex( 1,-1,-1),   // 3
                         GLVertex( 1, 1, 1),   // 4
                         GLVertex(-1, 1, 1),   // 5
                         GLVertex(-1,-1, 1),   // 6
                         GLVertex( 1,-1, 1)    // 7
                        };

    for ( int n=0;n<8;++n) {
    	root->vertex[n] = new GLVertex(*root->center + direction[n] * root->scale) ;
    	root->f[n] = -1;
    }
    root->bb.clear();
#ifdef MULTI_AXIS
// Multi Axis
    root->bb.addPoint(*root->center + GLVertex(-2.0,-2.0,-2.0) * root->scale); // caluclate the minimum x,y,z coordinates
    root->bb.addPoint(*root->center + GLVertex( 2.0, 2.0, 2.0) * root->scale); // caluclate the maximum x,y,z coordinates
#else
    root->bb.addPoint( *vertex[2] ); // vertex[2] has the minimum x,y,z coordinates
    root->bb.addPoint( *vertex[4] ); // vertex[4] has the max x,y,z
#endif
    root->childcount = 0;
}

void Octree::clearNode(Octnode* current) {
	if (((current->childcount) == 8))
		for(int m=0;m<8;++m)
			clearNode(current->child[m]);
	current->force_delete_children();
}

} // end namespace
// end of file octree.cpp
