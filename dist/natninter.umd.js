(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(factory((global.natninter = {})));
}(this, (function (exports) { 'use strict';

function createCommonjsModule(fn, module) {
	return module = { exports: {} }, fn(module, module.exports), module.exports;
}

var rhillVoronoiCore = createCommonjsModule(function (module) {
/*!
Copyright (C) 2010-2013 Raymond Hill: https://github.com/gorhill/Javascript-Voronoi
MIT License: See https://github.com/gorhill/Javascript-Voronoi/LICENSE.md
*/
/*
Author: Raymond Hill (rhill@raymondhill.net)
Contributor: Jesse Morgan (morgajel@gmail.com)
File: rhill-voronoi-core.js
Version: 0.98
Date: January 21, 2013
Description: This is my personal Javascript implementation of
Steven Fortune's algorithm to compute Voronoi diagrams.

License: See https://github.com/gorhill/Javascript-Voronoi/LICENSE.md
Credits: See https://github.com/gorhill/Javascript-Voronoi/CREDITS.md
History: See https://github.com/gorhill/Javascript-Voronoi/CHANGELOG.md

## Usage:

  var sites = [{x:300,y:300}, {x:100,y:100}, {x:200,y:500}, {x:250,y:450}, {x:600,y:150}];
  // xl, xr means x left, x right
  // yt, yb means y top, y bottom
  var bbox = {xl:0, xr:800, yt:0, yb:600};
  var voronoi = new Voronoi();
  // pass an object which exhibits xl, xr, yt, yb properties. The bounding
  // box will be used to connect unbound edges, and to close open cells
  result = voronoi.compute(sites, bbox);
  // render, further analyze, etc.

Return value:
  An object with the following properties:

  result.vertices = an array of unordered, unique Voronoi.Vertex objects making
    up the Voronoi diagram.
  result.edges = an array of unordered, unique Voronoi.Edge objects making up
    the Voronoi diagram.
  result.cells = an array of Voronoi.Cell object making up the Voronoi diagram.
    A Cell object might have an empty array of halfedges, meaning no Voronoi
    cell could be computed for a particular cell.
  result.execTime = the time it took to compute the Voronoi diagram, in
    milliseconds.

Voronoi.Vertex object:
  x: The x position of the vertex.
  y: The y position of the vertex.

Voronoi.Edge object:
  lSite: the Voronoi site object at the left of this Voronoi.Edge object.
  rSite: the Voronoi site object at the right of this Voronoi.Edge object (can
    be null).
  va: an object with an 'x' and a 'y' property defining the start point
    (relative to the Voronoi site on the left) of this Voronoi.Edge object.
  vb: an object with an 'x' and a 'y' property defining the end point
    (relative to Voronoi site on the left) of this Voronoi.Edge object.

  For edges which are used to close open cells (using the supplied bounding
  box), the rSite property will be null.

Voronoi.Cell object:
  site: the Voronoi site object associated with the Voronoi cell.
  halfedges: an array of Voronoi.Halfedge objects, ordered counterclockwise,
    defining the polygon for this Voronoi cell.

Voronoi.Halfedge object:
  site: the Voronoi site object owning this Voronoi.Halfedge object.
  edge: a reference to the unique Voronoi.Edge object underlying this
    Voronoi.Halfedge object.
  getStartpoint(): a method returning an object with an 'x' and a 'y' property
    for the start point of this halfedge. Keep in mind halfedges are always
    countercockwise.
  getEndpoint(): a method returning an object with an 'x' and a 'y' property
    for the end point of this halfedge. Keep in mind halfedges are always
    countercockwise.

TODO: Identify opportunities for performance improvement.

TODO: Let the user close the Voronoi cells, do not do it automatically. Not only let
      him close the cells, but also allow him to close more than once using a different
      bounding box for the same Voronoi diagram.
*/

/*global Math */

// ---------------------------------------------------------------------------

function Voronoi() {
    this.vertices = null;
    this.edges = null;
    this.cells = null;
    this.toRecycle = null;
    this.beachsectionJunkyard = [];
    this.circleEventJunkyard = [];
    this.vertexJunkyard = [];
    this.edgeJunkyard = [];
    this.cellJunkyard = [];
    }

// ---------------------------------------------------------------------------

Voronoi.prototype.reset = function() {
    if (!this.beachline) {
        this.beachline = new this.RBTree();
        }
    // Move leftover beachsections to the beachsection junkyard.
    if (this.beachline.root) {
        var beachsection = this.beachline.getFirst(this.beachline.root);
        while (beachsection) {
            this.beachsectionJunkyard.push(beachsection); // mark for reuse
            beachsection = beachsection.rbNext;
            }
        }
    this.beachline.root = null;
    if (!this.circleEvents) {
        this.circleEvents = new this.RBTree();
        }
    this.circleEvents.root = this.firstCircleEvent = null;
    this.vertices = [];
    this.edges = [];
    this.cells = [];
    };

Voronoi.prototype.sqrt = Math.sqrt;
Voronoi.prototype.abs = Math.abs;
Voronoi.prototype.ε = Voronoi.ε = 1e-9;
Voronoi.prototype.invε = Voronoi.invε = 1.0 / Voronoi.ε;
Voronoi.prototype.equalWithEpsilon = function(a,b){return this.abs(a-b)<1e-9;};
Voronoi.prototype.greaterThanWithEpsilon = function(a,b){return a-b>1e-9;};
Voronoi.prototype.greaterThanOrEqualWithEpsilon = function(a,b){return b-a<1e-9;};
Voronoi.prototype.lessThanWithEpsilon = function(a,b){return b-a>1e-9;};
Voronoi.prototype.lessThanOrEqualWithEpsilon = function(a,b){return a-b<1e-9;};

// ---------------------------------------------------------------------------
// Red-Black tree code (based on C version of "rbtree" by Franck Bui-Huu
// https://github.com/fbuihuu/libtree/blob/master/rb.c

Voronoi.prototype.RBTree = function() {
    this.root = null;
    };

Voronoi.prototype.RBTree.prototype.rbInsertSuccessor = function(node, successor) {
    var parent;
    if (node) {
        // >>> rhill 2011-05-27: Performance: cache previous/next nodes
        successor.rbPrevious = node;
        successor.rbNext = node.rbNext;
        if (node.rbNext) {
            node.rbNext.rbPrevious = successor;
            }
        node.rbNext = successor;
        // <<<
        if (node.rbRight) {
            // in-place expansion of node.rbRight.getFirst();
            node = node.rbRight;
            while (node.rbLeft) {node = node.rbLeft;}
            node.rbLeft = successor;
            }
        else {
            node.rbRight = successor;
            }
        parent = node;
        }
    // rhill 2011-06-07: if node is null, successor must be inserted
    // to the left-most part of the tree
    else if (this.root) {
        node = this.getFirst(this.root);
        // >>> Performance: cache previous/next nodes
        successor.rbPrevious = null;
        successor.rbNext = node;
        node.rbPrevious = successor;
        // <<<
        node.rbLeft = successor;
        parent = node;
        }
    else {
        // >>> Performance: cache previous/next nodes
        successor.rbPrevious = successor.rbNext = null;
        // <<<
        this.root = successor;
        parent = null;
        }
    successor.rbLeft = successor.rbRight = null;
    successor.rbParent = parent;
    successor.rbRed = true;
    // Fixup the modified tree by recoloring nodes and performing
    // rotations (2 at most) hence the red-black tree properties are
    // preserved.
    var grandpa, uncle;
    node = successor;
    while (parent && parent.rbRed) {
        grandpa = parent.rbParent;
        if (parent === grandpa.rbLeft) {
            uncle = grandpa.rbRight;
            if (uncle && uncle.rbRed) {
                parent.rbRed = uncle.rbRed = false;
                grandpa.rbRed = true;
                node = grandpa;
                }
            else {
                if (node === parent.rbRight) {
                    this.rbRotateLeft(parent);
                    node = parent;
                    parent = node.rbParent;
                    }
                parent.rbRed = false;
                grandpa.rbRed = true;
                this.rbRotateRight(grandpa);
                }
            }
        else {
            uncle = grandpa.rbLeft;
            if (uncle && uncle.rbRed) {
                parent.rbRed = uncle.rbRed = false;
                grandpa.rbRed = true;
                node = grandpa;
                }
            else {
                if (node === parent.rbLeft) {
                    this.rbRotateRight(parent);
                    node = parent;
                    parent = node.rbParent;
                    }
                parent.rbRed = false;
                grandpa.rbRed = true;
                this.rbRotateLeft(grandpa);
                }
            }
        parent = node.rbParent;
        }
    this.root.rbRed = false;
    };

Voronoi.prototype.RBTree.prototype.rbRemoveNode = function(node) {
    // >>> rhill 2011-05-27: Performance: cache previous/next nodes
    if (node.rbNext) {
        node.rbNext.rbPrevious = node.rbPrevious;
        }
    if (node.rbPrevious) {
        node.rbPrevious.rbNext = node.rbNext;
        }
    node.rbNext = node.rbPrevious = null;
    // <<<
    var parent = node.rbParent,
        left = node.rbLeft,
        right = node.rbRight,
        next;
    if (!left) {
        next = right;
        }
    else if (!right) {
        next = left;
        }
    else {
        next = this.getFirst(right);
        }
    if (parent) {
        if (parent.rbLeft === node) {
            parent.rbLeft = next;
            }
        else {
            parent.rbRight = next;
            }
        }
    else {
        this.root = next;
        }
    // enforce red-black rules
    var isRed;
    if (left && right) {
        isRed = next.rbRed;
        next.rbRed = node.rbRed;
        next.rbLeft = left;
        left.rbParent = next;
        if (next !== right) {
            parent = next.rbParent;
            next.rbParent = node.rbParent;
            node = next.rbRight;
            parent.rbLeft = node;
            next.rbRight = right;
            right.rbParent = next;
            }
        else {
            next.rbParent = parent;
            parent = next;
            node = next.rbRight;
            }
        }
    else {
        isRed = node.rbRed;
        node = next;
        }
    // 'node' is now the sole successor's child and 'parent' its
    // new parent (since the successor can have been moved)
    if (node) {
        node.rbParent = parent;
        }
    // the 'easy' cases
    if (isRed) {return;}
    if (node && node.rbRed) {
        node.rbRed = false;
        return;
        }
    // the other cases
    var sibling;
    do {
        if (node === this.root) {
            break;
            }
        if (node === parent.rbLeft) {
            sibling = parent.rbRight;
            if (sibling.rbRed) {
                sibling.rbRed = false;
                parent.rbRed = true;
                this.rbRotateLeft(parent);
                sibling = parent.rbRight;
                }
            if ((sibling.rbLeft && sibling.rbLeft.rbRed) || (sibling.rbRight && sibling.rbRight.rbRed)) {
                if (!sibling.rbRight || !sibling.rbRight.rbRed) {
                    sibling.rbLeft.rbRed = false;
                    sibling.rbRed = true;
                    this.rbRotateRight(sibling);
                    sibling = parent.rbRight;
                    }
                sibling.rbRed = parent.rbRed;
                parent.rbRed = sibling.rbRight.rbRed = false;
                this.rbRotateLeft(parent);
                node = this.root;
                break;
                }
            }
        else {
            sibling = parent.rbLeft;
            if (sibling.rbRed) {
                sibling.rbRed = false;
                parent.rbRed = true;
                this.rbRotateRight(parent);
                sibling = parent.rbLeft;
                }
            if ((sibling.rbLeft && sibling.rbLeft.rbRed) || (sibling.rbRight && sibling.rbRight.rbRed)) {
                if (!sibling.rbLeft || !sibling.rbLeft.rbRed) {
                    sibling.rbRight.rbRed = false;
                    sibling.rbRed = true;
                    this.rbRotateLeft(sibling);
                    sibling = parent.rbLeft;
                    }
                sibling.rbRed = parent.rbRed;
                parent.rbRed = sibling.rbLeft.rbRed = false;
                this.rbRotateRight(parent);
                node = this.root;
                break;
                }
            }
        sibling.rbRed = true;
        node = parent;
        parent = parent.rbParent;
    } while (!node.rbRed);
    if (node) {node.rbRed = false;}
    };

Voronoi.prototype.RBTree.prototype.rbRotateLeft = function(node) {
    var p = node,
        q = node.rbRight, // can't be null
        parent = p.rbParent;
    if (parent) {
        if (parent.rbLeft === p) {
            parent.rbLeft = q;
            }
        else {
            parent.rbRight = q;
            }
        }
    else {
        this.root = q;
        }
    q.rbParent = parent;
    p.rbParent = q;
    p.rbRight = q.rbLeft;
    if (p.rbRight) {
        p.rbRight.rbParent = p;
        }
    q.rbLeft = p;
    };

Voronoi.prototype.RBTree.prototype.rbRotateRight = function(node) {
    var p = node,
        q = node.rbLeft, // can't be null
        parent = p.rbParent;
    if (parent) {
        if (parent.rbLeft === p) {
            parent.rbLeft = q;
            }
        else {
            parent.rbRight = q;
            }
        }
    else {
        this.root = q;
        }
    q.rbParent = parent;
    p.rbParent = q;
    p.rbLeft = q.rbRight;
    if (p.rbLeft) {
        p.rbLeft.rbParent = p;
        }
    q.rbRight = p;
    };

Voronoi.prototype.RBTree.prototype.getFirst = function(node) {
    while (node.rbLeft) {
        node = node.rbLeft;
        }
    return node;
    };

Voronoi.prototype.RBTree.prototype.getLast = function(node) {
    while (node.rbRight) {
        node = node.rbRight;
        }
    return node;
    };

// ---------------------------------------------------------------------------
// Diagram methods

Voronoi.prototype.Diagram = function(site) {
    this.site = site;
    };

// ---------------------------------------------------------------------------
// Cell methods

Voronoi.prototype.Cell = function(site) {
    this.site = site;
    this.halfedges = [];
    this.closeMe = false;
    };

Voronoi.prototype.Cell.prototype.init = function(site) {
    this.site = site;
    this.halfedges = [];
    this.closeMe = false;
    return this;
    };

Voronoi.prototype.createCell = function(site) {
    var cell = this.cellJunkyard.pop();
    if ( cell ) {
        return cell.init(site);
        }
    return new this.Cell(site);
    };

Voronoi.prototype.Cell.prototype.prepareHalfedges = function() {
    var halfedges = this.halfedges,
        iHalfedge = halfedges.length,
        edge;
    // get rid of unused halfedges
    // rhill 2011-05-27: Keep it simple, no point here in trying
    // to be fancy: dangling edges are a typically a minority.
    while (iHalfedge--) {
        edge = halfedges[iHalfedge].edge;
        if (!edge.vb || !edge.va) {
            halfedges.splice(iHalfedge,1);
            }
        }

    // rhill 2011-05-26: I tried to use a binary search at insertion
    // time to keep the array sorted on-the-fly (in Cell.addHalfedge()).
    // There was no real benefits in doing so, performance on
    // Firefox 3.6 was improved marginally, while performance on
    // Opera 11 was penalized marginally.
    halfedges.sort(function(a,b){return b.angle-a.angle;});
    return halfedges.length;
    };

// Return a list of the neighbor Ids
Voronoi.prototype.Cell.prototype.getNeighborIds = function() {
    var neighbors = [],
        iHalfedge = this.halfedges.length,
        edge;
    while (iHalfedge--){
        edge = this.halfedges[iHalfedge].edge;
        if (edge.lSite !== null && edge.lSite.voronoiId != this.site.voronoiId) {
            neighbors.push(edge.lSite.voronoiId);
            }
        else if (edge.rSite !== null && edge.rSite.voronoiId != this.site.voronoiId){
            neighbors.push(edge.rSite.voronoiId);
            }
        }
    return neighbors;
    };

// Compute bounding box
//
Voronoi.prototype.Cell.prototype.getBbox = function() {
    var halfedges = this.halfedges,
        iHalfedge = halfedges.length,
        xmin = Infinity,
        ymin = Infinity,
        xmax = -Infinity,
        ymax = -Infinity,
        v, vx, vy;
    while (iHalfedge--) {
        v = halfedges[iHalfedge].getStartpoint();
        vx = v.x;
        vy = v.y;
        if (vx < xmin) {xmin = vx;}
        if (vy < ymin) {ymin = vy;}
        if (vx > xmax) {xmax = vx;}
        if (vy > ymax) {ymax = vy;}
        // we dont need to take into account end point,
        // since each end point matches a start point
        }
    return {
        x: xmin,
        y: ymin,
        width: xmax-xmin,
        height: ymax-ymin
        };
    };

// Return whether a point is inside, on, or outside the cell:
//   -1: point is outside the perimeter of the cell
//    0: point is on the perimeter of the cell
//    1: point is inside the perimeter of the cell
//
Voronoi.prototype.Cell.prototype.pointIntersection = function(x, y) {
    // Check if point in polygon. Since all polygons of a Voronoi
    // diagram are convex, then:
    // http://paulbourke.net/geometry/polygonmesh/
    // Solution 3 (2D):
    //   "If the polygon is convex then one can consider the polygon
    //   "as a 'path' from the first vertex. A point is on the interior
    //   "of this polygons if it is always on the same side of all the
    //   "line segments making up the path. ...
    //   "(y - y0) (x1 - x0) - (x - x0) (y1 - y0)
    //   "if it is less than 0 then P is to the right of the line segment,
    //   "if greater than 0 it is to the left, if equal to 0 then it lies
    //   "on the line segment"
    var halfedges = this.halfedges,
        iHalfedge = halfedges.length,
        halfedge,
        p0, p1, r;
    while (iHalfedge--) {
        halfedge = halfedges[iHalfedge];
        p0 = halfedge.getStartpoint();
        p1 = halfedge.getEndpoint();
        r = (y-p0.y)*(p1.x-p0.x)-(x-p0.x)*(p1.y-p0.y);
        if (!r) {
            return 0;
            }
        if (r > 0) {
            return -1;
            }
        }
    return 1;
    };

// ---------------------------------------------------------------------------
// Edge methods
//

Voronoi.prototype.Vertex = function(x, y) {
    this.x = x;
    this.y = y;
    };

Voronoi.prototype.Edge = function(lSite, rSite) {
    this.lSite = lSite;
    this.rSite = rSite;
    this.va = this.vb = null;
    };

Voronoi.prototype.Halfedge = function(edge, lSite, rSite) {
    this.site = lSite;
    this.edge = edge;
    // 'angle' is a value to be used for properly sorting the
    // halfsegments counterclockwise. By convention, we will
    // use the angle of the line defined by the 'site to the left'
    // to the 'site to the right'.
    // However, border edges have no 'site to the right': thus we
    // use the angle of line perpendicular to the halfsegment (the
    // edge should have both end points defined in such case.)
    if (rSite) {
        this.angle = Math.atan2(rSite.y-lSite.y, rSite.x-lSite.x);
        }
    else {
        var va = edge.va,
            vb = edge.vb;
        // rhill 2011-05-31: used to call getStartpoint()/getEndpoint(),
        // but for performance purpose, these are expanded in place here.
        this.angle = edge.lSite === lSite ?
            Math.atan2(vb.x-va.x, va.y-vb.y) :
            Math.atan2(va.x-vb.x, vb.y-va.y);
        }
    };

Voronoi.prototype.createHalfedge = function(edge, lSite, rSite) {
    return new this.Halfedge(edge, lSite, rSite);
    };

Voronoi.prototype.Halfedge.prototype.getStartpoint = function() {
    return this.edge.lSite === this.site ? this.edge.va : this.edge.vb;
    };

Voronoi.prototype.Halfedge.prototype.getEndpoint = function() {
    return this.edge.lSite === this.site ? this.edge.vb : this.edge.va;
    };



// this create and add a vertex to the internal collection

Voronoi.prototype.createVertex = function(x, y) {
    var v = this.vertexJunkyard.pop();
    if ( !v ) {
        v = new this.Vertex(x, y);
        }
    else {
        v.x = x;
        v.y = y;
        }
    this.vertices.push(v);
    return v;
    };

// this create and add an edge to internal collection, and also create
// two halfedges which are added to each site's counterclockwise array
// of halfedges.

Voronoi.prototype.createEdge = function(lSite, rSite, va, vb) {
    var edge = this.edgeJunkyard.pop();
    if ( !edge ) {
        edge = new this.Edge(lSite, rSite);
        }
    else {
        edge.lSite = lSite;
        edge.rSite = rSite;
        edge.va = edge.vb = null;
        }

    this.edges.push(edge);
    if (va) {
        this.setEdgeStartpoint(edge, lSite, rSite, va);
        }
    if (vb) {
        this.setEdgeEndpoint(edge, lSite, rSite, vb);
        }
    this.cells[lSite.voronoiId].halfedges.push(this.createHalfedge(edge, lSite, rSite));
    this.cells[rSite.voronoiId].halfedges.push(this.createHalfedge(edge, rSite, lSite));
    return edge;
    };

Voronoi.prototype.createBorderEdge = function(lSite, va, vb) {
    var edge = this.edgeJunkyard.pop();
    if ( !edge ) {
        edge = new this.Edge(lSite, null);
        }
    else {
        edge.lSite = lSite;
        edge.rSite = null;
        }
    edge.va = va;
    edge.vb = vb;
    this.edges.push(edge);
    return edge;
    };

Voronoi.prototype.setEdgeStartpoint = function(edge, lSite, rSite, vertex) {
    if (!edge.va && !edge.vb) {
        edge.va = vertex;
        edge.lSite = lSite;
        edge.rSite = rSite;
        }
    else if (edge.lSite === rSite) {
        edge.vb = vertex;
        }
    else {
        edge.va = vertex;
        }
    };

Voronoi.prototype.setEdgeEndpoint = function(edge, lSite, rSite, vertex) {
    this.setEdgeStartpoint(edge, rSite, lSite, vertex);
    };

// ---------------------------------------------------------------------------
// Beachline methods

// rhill 2011-06-07: For some reasons, performance suffers significantly
// when instanciating a literal object instead of an empty ctor
Voronoi.prototype.Beachsection = function() {
    };

// rhill 2011-06-02: A lot of Beachsection instanciations
// occur during the computation of the Voronoi diagram,
// somewhere between the number of sites and twice the
// number of sites, while the number of Beachsections on the
// beachline at any given time is comparatively low. For this
// reason, we reuse already created Beachsections, in order
// to avoid new memory allocation. This resulted in a measurable
// performance gain.

Voronoi.prototype.createBeachsection = function(site) {
    var beachsection = this.beachsectionJunkyard.pop();
    if (!beachsection) {
        beachsection = new this.Beachsection();
        }
    beachsection.site = site;
    return beachsection;
    };

// calculate the left break point of a particular beach section,
// given a particular sweep line
Voronoi.prototype.leftBreakPoint = function(arc, directrix) {
    // http://en.wikipedia.org/wiki/Parabola
    // http://en.wikipedia.org/wiki/Quadratic_equation
    // h1 = x1,
    // k1 = (y1+directrix)/2,
    // h2 = x2,
    // k2 = (y2+directrix)/2,
    // p1 = k1-directrix,
    // a1 = 1/(4*p1),
    // b1 = -h1/(2*p1),
    // c1 = h1*h1/(4*p1)+k1,
    // p2 = k2-directrix,
    // a2 = 1/(4*p2),
    // b2 = -h2/(2*p2),
    // c2 = h2*h2/(4*p2)+k2,
    // x = (-(b2-b1) + Math.sqrt((b2-b1)*(b2-b1) - 4*(a2-a1)*(c2-c1))) / (2*(a2-a1))
    // When x1 become the x-origin:
    // h1 = 0,
    // k1 = (y1+directrix)/2,
    // h2 = x2-x1,
    // k2 = (y2+directrix)/2,
    // p1 = k1-directrix,
    // a1 = 1/(4*p1),
    // b1 = 0,
    // c1 = k1,
    // p2 = k2-directrix,
    // a2 = 1/(4*p2),
    // b2 = -h2/(2*p2),
    // c2 = h2*h2/(4*p2)+k2,
    // x = (-b2 + Math.sqrt(b2*b2 - 4*(a2-a1)*(c2-k1))) / (2*(a2-a1)) + x1

    // change code below at your own risk: care has been taken to
    // reduce errors due to computers' finite arithmetic precision.
    // Maybe can still be improved, will see if any more of this
    // kind of errors pop up again.
    var site = arc.site,
        rfocx = site.x,
        rfocy = site.y,
        pby2 = rfocy-directrix;
    // parabola in degenerate case where focus is on directrix
    if (!pby2) {
        return rfocx;
        }
    var lArc = arc.rbPrevious;
    if (!lArc) {
        return -Infinity;
        }
    site = lArc.site;
    var lfocx = site.x,
        lfocy = site.y,
        plby2 = lfocy-directrix;
    // parabola in degenerate case where focus is on directrix
    if (!plby2) {
        return lfocx;
        }
    var hl = lfocx-rfocx,
        aby2 = 1/pby2-1/plby2,
        b = hl/plby2;
    if (aby2) {
        return (-b+this.sqrt(b*b-2*aby2*(hl*hl/(-2*plby2)-lfocy+plby2/2+rfocy-pby2/2)))/aby2+rfocx;
        }
    // both parabolas have same distance to directrix, thus break point is midway
    return (rfocx+lfocx)/2;
    };

// calculate the right break point of a particular beach section,
// given a particular directrix
Voronoi.prototype.rightBreakPoint = function(arc, directrix) {
    var rArc = arc.rbNext;
    if (rArc) {
        return this.leftBreakPoint(rArc, directrix);
        }
    var site = arc.site;
    return site.y === directrix ? site.x : Infinity;
    };

Voronoi.prototype.detachBeachsection = function(beachsection) {
    this.detachCircleEvent(beachsection); // detach potentially attached circle event
    this.beachline.rbRemoveNode(beachsection); // remove from RB-tree
    this.beachsectionJunkyard.push(beachsection); // mark for reuse
    };

Voronoi.prototype.removeBeachsection = function(beachsection) {
    var circle = beachsection.circleEvent,
        x = circle.x,
        y = circle.ycenter,
        vertex = this.createVertex(x, y),
        previous = beachsection.rbPrevious,
        next = beachsection.rbNext,
        disappearingTransitions = [beachsection],
        abs_fn = Math.abs;

    // remove collapsed beachsection from beachline
    this.detachBeachsection(beachsection);

    // there could be more than one empty arc at the deletion point, this
    // happens when more than two edges are linked by the same vertex,
    // so we will collect all those edges by looking up both sides of
    // the deletion point.
    // by the way, there is *always* a predecessor/successor to any collapsed
    // beach section, it's just impossible to have a collapsing first/last
    // beach sections on the beachline, since they obviously are unconstrained
    // on their left/right side.

    // look left
    var lArc = previous;
    while (lArc.circleEvent && abs_fn(x-lArc.circleEvent.x)<1e-9 && abs_fn(y-lArc.circleEvent.ycenter)<1e-9) {
        previous = lArc.rbPrevious;
        disappearingTransitions.unshift(lArc);
        this.detachBeachsection(lArc); // mark for reuse
        lArc = previous;
        }
    // even though it is not disappearing, I will also add the beach section
    // immediately to the left of the left-most collapsed beach section, for
    // convenience, since we need to refer to it later as this beach section
    // is the 'left' site of an edge for which a start point is set.
    disappearingTransitions.unshift(lArc);
    this.detachCircleEvent(lArc);

    // look right
    var rArc = next;
    while (rArc.circleEvent && abs_fn(x-rArc.circleEvent.x)<1e-9 && abs_fn(y-rArc.circleEvent.ycenter)<1e-9) {
        next = rArc.rbNext;
        disappearingTransitions.push(rArc);
        this.detachBeachsection(rArc); // mark for reuse
        rArc = next;
        }
    // we also have to add the beach section immediately to the right of the
    // right-most collapsed beach section, since there is also a disappearing
    // transition representing an edge's start point on its left.
    disappearingTransitions.push(rArc);
    this.detachCircleEvent(rArc);

    // walk through all the disappearing transitions between beach sections and
    // set the start point of their (implied) edge.
    var nArcs = disappearingTransitions.length,
        iArc;
    for (iArc=1; iArc<nArcs; iArc++) {
        rArc = disappearingTransitions[iArc];
        lArc = disappearingTransitions[iArc-1];
        this.setEdgeStartpoint(rArc.edge, lArc.site, rArc.site, vertex);
        }

    // create a new edge as we have now a new transition between
    // two beach sections which were previously not adjacent.
    // since this edge appears as a new vertex is defined, the vertex
    // actually define an end point of the edge (relative to the site
    // on the left)
    lArc = disappearingTransitions[0];
    rArc = disappearingTransitions[nArcs-1];
    rArc.edge = this.createEdge(lArc.site, rArc.site, undefined, vertex);

    // create circle events if any for beach sections left in the beachline
    // adjacent to collapsed sections
    this.attachCircleEvent(lArc);
    this.attachCircleEvent(rArc);
    };

Voronoi.prototype.addBeachsection = function(site) {
    var x = site.x,
        directrix = site.y;

    // find the left and right beach sections which will surround the newly
    // created beach section.
    // rhill 2011-06-01: This loop is one of the most often executed,
    // hence we expand in-place the comparison-against-epsilon calls.
    var lArc, rArc,
        dxl, dxr,
        node = this.beachline.root;

    while (node) {
        dxl = this.leftBreakPoint(node,directrix)-x;
        // x lessThanWithEpsilon xl => falls somewhere before the left edge of the beachsection
        if (dxl > 1e-9) {
            // this case should never happen
            // if (!node.rbLeft) {
            //    rArc = node.rbLeft;
            //    break;
            //    }
            node = node.rbLeft;
            }
        else {
            dxr = x-this.rightBreakPoint(node,directrix);
            // x greaterThanWithEpsilon xr => falls somewhere after the right edge of the beachsection
            if (dxr > 1e-9) {
                if (!node.rbRight) {
                    lArc = node;
                    break;
                    }
                node = node.rbRight;
                }
            else {
                // x equalWithEpsilon xl => falls exactly on the left edge of the beachsection
                if (dxl > -1e-9) {
                    lArc = node.rbPrevious;
                    rArc = node;
                    }
                // x equalWithEpsilon xr => falls exactly on the right edge of the beachsection
                else if (dxr > -1e-9) {
                    lArc = node;
                    rArc = node.rbNext;
                    }
                // falls exactly somewhere in the middle of the beachsection
                else {
                    lArc = rArc = node;
                    }
                break;
                }
            }
        }
    // at this point, keep in mind that lArc and/or rArc could be
    // undefined or null.

    // create a new beach section object for the site and add it to RB-tree
    var newArc = this.createBeachsection(site);
    this.beachline.rbInsertSuccessor(lArc, newArc);

    // cases:
    //

    // [null,null]
    // least likely case: new beach section is the first beach section on the
    // beachline.
    // This case means:
    //   no new transition appears
    //   no collapsing beach section
    //   new beachsection become root of the RB-tree
    if (!lArc && !rArc) {
        return;
        }

    // [lArc,rArc] where lArc == rArc
    // most likely case: new beach section split an existing beach
    // section.
    // This case means:
    //   one new transition appears
    //   the left and right beach section might be collapsing as a result
    //   two new nodes added to the RB-tree
    if (lArc === rArc) {
        // invalidate circle event of split beach section
        this.detachCircleEvent(lArc);

        // split the beach section into two separate beach sections
        rArc = this.createBeachsection(lArc.site);
        this.beachline.rbInsertSuccessor(newArc, rArc);

        // since we have a new transition between two beach sections,
        // a new edge is born
        newArc.edge = rArc.edge = this.createEdge(lArc.site, newArc.site);

        // check whether the left and right beach sections are collapsing
        // and if so create circle events, to be notified when the point of
        // collapse is reached.
        this.attachCircleEvent(lArc);
        this.attachCircleEvent(rArc);
        return;
        }

    // [lArc,null]
    // even less likely case: new beach section is the *last* beach section
    // on the beachline -- this can happen *only* if *all* the previous beach
    // sections currently on the beachline share the same y value as
    // the new beach section.
    // This case means:
    //   one new transition appears
    //   no collapsing beach section as a result
    //   new beach section become right-most node of the RB-tree
    if (lArc && !rArc) {
        newArc.edge = this.createEdge(lArc.site,newArc.site);
        return;
        }

    // [null,rArc]
    // impossible case: because sites are strictly processed from top to bottom,
    // and left to right, which guarantees that there will always be a beach section
    // on the left -- except of course when there are no beach section at all on
    // the beach line, which case was handled above.
    // rhill 2011-06-02: No point testing in non-debug version
    //if (!lArc && rArc) {
    //    throw "Voronoi.addBeachsection(): What is this I don't even";
    //    }

    // [lArc,rArc] where lArc != rArc
    // somewhat less likely case: new beach section falls *exactly* in between two
    // existing beach sections
    // This case means:
    //   one transition disappears
    //   two new transitions appear
    //   the left and right beach section might be collapsing as a result
    //   only one new node added to the RB-tree
    if (lArc !== rArc) {
        // invalidate circle events of left and right sites
        this.detachCircleEvent(lArc);
        this.detachCircleEvent(rArc);

        // an existing transition disappears, meaning a vertex is defined at
        // the disappearance point.
        // since the disappearance is caused by the new beachsection, the
        // vertex is at the center of the circumscribed circle of the left,
        // new and right beachsections.
        // http://mathforum.org/library/drmath/view/55002.html
        // Except that I bring the origin at A to simplify
        // calculation
        var lSite = lArc.site,
            ax = lSite.x,
            ay = lSite.y,
            bx=site.x-ax,
            by=site.y-ay,
            rSite = rArc.site,
            cx=rSite.x-ax,
            cy=rSite.y-ay,
            d=2*(bx*cy-by*cx),
            hb=bx*bx+by*by,
            hc=cx*cx+cy*cy,
            vertex = this.createVertex((cy*hb-by*hc)/d+ax, (bx*hc-cx*hb)/d+ay);

        // one transition disappear
        this.setEdgeStartpoint(rArc.edge, lSite, rSite, vertex);

        // two new transitions appear at the new vertex location
        newArc.edge = this.createEdge(lSite, site, undefined, vertex);
        rArc.edge = this.createEdge(site, rSite, undefined, vertex);

        // check whether the left and right beach sections are collapsing
        // and if so create circle events, to handle the point of collapse.
        this.attachCircleEvent(lArc);
        this.attachCircleEvent(rArc);
        return;
        }
    };

// ---------------------------------------------------------------------------
// Circle event methods

// rhill 2011-06-07: For some reasons, performance suffers significantly
// when instanciating a literal object instead of an empty ctor
Voronoi.prototype.CircleEvent = function() {
    // rhill 2013-10-12: it helps to state exactly what we are at ctor time.
    this.arc = null;
    this.rbLeft = null;
    this.rbNext = null;
    this.rbParent = null;
    this.rbPrevious = null;
    this.rbRed = false;
    this.rbRight = null;
    this.site = null;
    this.x = this.y = this.ycenter = 0;
    };

Voronoi.prototype.attachCircleEvent = function(arc) {
    var lArc = arc.rbPrevious,
        rArc = arc.rbNext;
    if (!lArc || !rArc) {return;} // does that ever happen?
    var lSite = lArc.site,
        cSite = arc.site,
        rSite = rArc.site;

    // If site of left beachsection is same as site of
    // right beachsection, there can't be convergence
    if (lSite===rSite) {return;}

    // Find the circumscribed circle for the three sites associated
    // with the beachsection triplet.
    // rhill 2011-05-26: It is more efficient to calculate in-place
    // rather than getting the resulting circumscribed circle from an
    // object returned by calling Voronoi.circumcircle()
    // http://mathforum.org/library/drmath/view/55002.html
    // Except that I bring the origin at cSite to simplify calculations.
    // The bottom-most part of the circumcircle is our Fortune 'circle
    // event', and its center is a vertex potentially part of the final
    // Voronoi diagram.
    var bx = cSite.x,
        by = cSite.y,
        ax = lSite.x-bx,
        ay = lSite.y-by,
        cx = rSite.x-bx,
        cy = rSite.y-by;

    // If points l->c->r are clockwise, then center beach section does not
    // collapse, hence it can't end up as a vertex (we reuse 'd' here, which
    // sign is reverse of the orientation, hence we reverse the test.
    // http://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon
    // rhill 2011-05-21: Nasty finite precision error which caused circumcircle() to
    // return infinites: 1e-12 seems to fix the problem.
    var d = 2*(ax*cy-ay*cx);
    if (d >= -2e-12){return;}

    var ha = ax*ax+ay*ay,
        hc = cx*cx+cy*cy,
        x = (cy*ha-ay*hc)/d,
        y = (ax*hc-cx*ha)/d,
        ycenter = y+by;

    // Important: ybottom should always be under or at sweep, so no need
    // to waste CPU cycles by checking

    // recycle circle event object if possible
    var circleEvent = this.circleEventJunkyard.pop();
    if (!circleEvent) {
        circleEvent = new this.CircleEvent();
        }
    circleEvent.arc = arc;
    circleEvent.site = cSite;
    circleEvent.x = x+bx;
    circleEvent.y = ycenter+this.sqrt(x*x+y*y); // y bottom
    circleEvent.ycenter = ycenter;
    arc.circleEvent = circleEvent;

    // find insertion point in RB-tree: circle events are ordered from
    // smallest to largest
    var predecessor = null,
        node = this.circleEvents.root;
    while (node) {
        if (circleEvent.y < node.y || (circleEvent.y === node.y && circleEvent.x <= node.x)) {
            if (node.rbLeft) {
                node = node.rbLeft;
                }
            else {
                predecessor = node.rbPrevious;
                break;
                }
            }
        else {
            if (node.rbRight) {
                node = node.rbRight;
                }
            else {
                predecessor = node;
                break;
                }
            }
        }
    this.circleEvents.rbInsertSuccessor(predecessor, circleEvent);
    if (!predecessor) {
        this.firstCircleEvent = circleEvent;
        }
    };

Voronoi.prototype.detachCircleEvent = function(arc) {
    var circleEvent = arc.circleEvent;
    if (circleEvent) {
        if (!circleEvent.rbPrevious) {
            this.firstCircleEvent = circleEvent.rbNext;
            }
        this.circleEvents.rbRemoveNode(circleEvent); // remove from RB-tree
        this.circleEventJunkyard.push(circleEvent);
        arc.circleEvent = null;
        }
    };

// ---------------------------------------------------------------------------
// Diagram completion methods

// connect dangling edges (not if a cursory test tells us
// it is not going to be visible.
// return value:
//   false: the dangling endpoint couldn't be connected
//   true: the dangling endpoint could be connected
Voronoi.prototype.connectEdge = function(edge, bbox) {
    // skip if end point already connected
    var vb = edge.vb;
    if (!!vb) {return true;}

    // make local copy for performance purpose
    var va = edge.va,
        xl = bbox.xl,
        xr = bbox.xr,
        yt = bbox.yt,
        yb = bbox.yb,
        lSite = edge.lSite,
        rSite = edge.rSite,
        lx = lSite.x,
        ly = lSite.y,
        rx = rSite.x,
        ry = rSite.y,
        fx = (lx+rx)/2,
        fy = (ly+ry)/2,
        fm, fb;

    // if we reach here, this means cells which use this edge will need
    // to be closed, whether because the edge was removed, or because it
    // was connected to the bounding box.
    this.cells[lSite.voronoiId].closeMe = true;
    this.cells[rSite.voronoiId].closeMe = true;

    // get the line equation of the bisector if line is not vertical
    if (ry !== ly) {
        fm = (lx-rx)/(ry-ly);
        fb = fy-fm*fx;
        }

    // remember, direction of line (relative to left site):
    // upward: left.x < right.x
    // downward: left.x > right.x
    // horizontal: left.x == right.x
    // upward: left.x < right.x
    // rightward: left.y < right.y
    // leftward: left.y > right.y
    // vertical: left.y == right.y

    // depending on the direction, find the best side of the
    // bounding box to use to determine a reasonable start point

    // rhill 2013-12-02:
    // While at it, since we have the values which define the line,
    // clip the end of va if it is outside the bbox.
    // https://github.com/gorhill/Javascript-Voronoi/issues/15
    // TODO: Do all the clipping here rather than rely on Liang-Barsky
    // which does not do well sometimes due to loss of arithmetic
    // precision. The code here doesn't degrade if one of the vertex is
    // at a huge distance.

    // special case: vertical line
    if (fm === undefined) {
        // doesn't intersect with viewport
        if (fx < xl || fx >= xr) {return false;}
        // downward
        if (lx > rx) {
            if (!va || va.y < yt) {
                va = this.createVertex(fx, yt);
                }
            else if (va.y >= yb) {
                return false;
                }
            vb = this.createVertex(fx, yb);
            }
        // upward
        else {
            if (!va || va.y > yb) {
                va = this.createVertex(fx, yb);
                }
            else if (va.y < yt) {
                return false;
                }
            vb = this.createVertex(fx, yt);
            }
        }
    // closer to vertical than horizontal, connect start point to the
    // top or bottom side of the bounding box
    else if (fm < -1 || fm > 1) {
        // downward
        if (lx > rx) {
            if (!va || va.y < yt) {
                va = this.createVertex((yt-fb)/fm, yt);
                }
            else if (va.y >= yb) {
                return false;
                }
            vb = this.createVertex((yb-fb)/fm, yb);
            }
        // upward
        else {
            if (!va || va.y > yb) {
                va = this.createVertex((yb-fb)/fm, yb);
                }
            else if (va.y < yt) {
                return false;
                }
            vb = this.createVertex((yt-fb)/fm, yt);
            }
        }
    // closer to horizontal than vertical, connect start point to the
    // left or right side of the bounding box
    else {
        // rightward
        if (ly < ry) {
            if (!va || va.x < xl) {
                va = this.createVertex(xl, fm*xl+fb);
                }
            else if (va.x >= xr) {
                return false;
                }
            vb = this.createVertex(xr, fm*xr+fb);
            }
        // leftward
        else {
            if (!va || va.x > xr) {
                va = this.createVertex(xr, fm*xr+fb);
                }
            else if (va.x < xl) {
                return false;
                }
            vb = this.createVertex(xl, fm*xl+fb);
            }
        }
    edge.va = va;
    edge.vb = vb;

    return true;
    };

// line-clipping code taken from:
//   Liang-Barsky function by Daniel White
//   http://www.skytopia.com/project/articles/compsci/clipping.html
// Thanks!
// A bit modified to minimize code paths
Voronoi.prototype.clipEdge = function(edge, bbox) {
    var ax = edge.va.x,
        ay = edge.va.y,
        bx = edge.vb.x,
        by = edge.vb.y,
        t0 = 0,
        t1 = 1,
        dx = bx-ax,
        dy = by-ay;
    // left
    var q = ax-bbox.xl;
    if (dx===0 && q<0) {return false;}
    var r = -q/dx;
    if (dx<0) {
        if (r<t0) {return false;}
        if (r<t1) {t1=r;}
        }
    else if (dx>0) {
        if (r>t1) {return false;}
        if (r>t0) {t0=r;}
        }
    // right
    q = bbox.xr-ax;
    if (dx===0 && q<0) {return false;}
    r = q/dx;
    if (dx<0) {
        if (r>t1) {return false;}
        if (r>t0) {t0=r;}
        }
    else if (dx>0) {
        if (r<t0) {return false;}
        if (r<t1) {t1=r;}
        }
    // top
    q = ay-bbox.yt;
    if (dy===0 && q<0) {return false;}
    r = -q/dy;
    if (dy<0) {
        if (r<t0) {return false;}
        if (r<t1) {t1=r;}
        }
    else if (dy>0) {
        if (r>t1) {return false;}
        if (r>t0) {t0=r;}
        }
    // bottom        
    q = bbox.yb-ay;
    if (dy===0 && q<0) {return false;}
    r = q/dy;
    if (dy<0) {
        if (r>t1) {return false;}
        if (r>t0) {t0=r;}
        }
    else if (dy>0) {
        if (r<t0) {return false;}
        if (r<t1) {t1=r;}
        }

    // if we reach this point, Voronoi edge is within bbox

    // if t0 > 0, va needs to change
    // rhill 2011-06-03: we need to create a new vertex rather
    // than modifying the existing one, since the existing
    // one is likely shared with at least another edge
    if (t0 > 0) {
        edge.va = this.createVertex(ax+t0*dx, ay+t0*dy);
        }

    // if t1 < 1, vb needs to change
    // rhill 2011-06-03: we need to create a new vertex rather
    // than modifying the existing one, since the existing
    // one is likely shared with at least another edge
    if (t1 < 1) {
        edge.vb = this.createVertex(ax+t1*dx, ay+t1*dy);
        }

    // va and/or vb were clipped, thus we will need to close
    // cells which use this edge.
    if ( t0 > 0 || t1 < 1 ) {
        this.cells[edge.lSite.voronoiId].closeMe = true;
        this.cells[edge.rSite.voronoiId].closeMe = true;
    }

    return true;
    };

// Connect/cut edges at bounding box
Voronoi.prototype.clipEdges = function(bbox) {
    // connect all dangling edges to bounding box
    // or get rid of them if it can't be done
    var edges = this.edges,
        iEdge = edges.length,
        edge,
        abs_fn = Math.abs;

    // iterate backward so we can splice safely
    while (iEdge--) {
        edge = edges[iEdge];
        // edge is removed if:
        //   it is wholly outside the bounding box
        //   it is looking more like a point than a line
        if (!this.connectEdge(edge, bbox) ||
            !this.clipEdge(edge, bbox) ||
            (abs_fn(edge.va.x-edge.vb.x)<1e-9 && abs_fn(edge.va.y-edge.vb.y)<1e-9)) {
            edge.va = edge.vb = null;
            edges.splice(iEdge,1);
            }
        }
    };

// Close the cells.
// The cells are bound by the supplied bounding box.
// Each cell refers to its associated site, and a list
// of halfedges ordered counterclockwise.
Voronoi.prototype.closeCells = function(bbox) {
    var xl = bbox.xl,
        xr = bbox.xr,
        yt = bbox.yt,
        yb = bbox.yb,
        cells = this.cells,
        iCell = cells.length,
        cell,
        iLeft,
        halfedges, nHalfedges,
        edge,
        va, vb, vz,
        lastBorderSegment,
        abs_fn = Math.abs;

    while (iCell--) {
        cell = cells[iCell];
        // prune, order halfedges counterclockwise, then add missing ones
        // required to close cells
        if (!cell.prepareHalfedges()) {
            continue;
            }
        if (!cell.closeMe) {
            continue;
            }
        // find first 'unclosed' point.
        // an 'unclosed' point will be the end point of a halfedge which
        // does not match the start point of the following halfedge
        halfedges = cell.halfedges;
        nHalfedges = halfedges.length;
        // special case: only one site, in which case, the viewport is the cell
        // ...

        // all other cases
        iLeft = 0;
        while (iLeft < nHalfedges) {
            va = halfedges[iLeft].getEndpoint();
            vz = halfedges[(iLeft+1) % nHalfedges].getStartpoint();
            // if end point is not equal to start point, we need to add the missing
            // halfedge(s) up to vz
            if (abs_fn(va.x-vz.x)>=1e-9 || abs_fn(va.y-vz.y)>=1e-9) {

                // rhill 2013-12-02:
                // "Holes" in the halfedges are not necessarily always adjacent.
                // https://github.com/gorhill/Javascript-Voronoi/issues/16

                // find entry point:
                switch (true) {

                    // walk downward along left side
                    case this.equalWithEpsilon(va.x,xl) && this.lessThanWithEpsilon(va.y,yb):
                        lastBorderSegment = this.equalWithEpsilon(vz.x,xl);
                        vb = this.createVertex(xl, lastBorderSegment ? vz.y : yb);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                    // walk rightward along bottom side
                    case this.equalWithEpsilon(va.y,yb) && this.lessThanWithEpsilon(va.x,xr):
                        lastBorderSegment = this.equalWithEpsilon(vz.y,yb);
                        vb = this.createVertex(lastBorderSegment ? vz.x : xr, yb);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                    // walk upward along right side
                    case this.equalWithEpsilon(va.x,xr) && this.greaterThanWithEpsilon(va.y,yt):
                        lastBorderSegment = this.equalWithEpsilon(vz.x,xr);
                        vb = this.createVertex(xr, lastBorderSegment ? vz.y : yt);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                    // walk leftward along top side
                    case this.equalWithEpsilon(va.y,yt) && this.greaterThanWithEpsilon(va.x,xl):
                        lastBorderSegment = this.equalWithEpsilon(vz.y,yt);
                        vb = this.createVertex(lastBorderSegment ? vz.x : xl, yt);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                        // walk downward along left side
                        lastBorderSegment = this.equalWithEpsilon(vz.x,xl);
                        vb = this.createVertex(xl, lastBorderSegment ? vz.y : yb);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                        // walk rightward along bottom side
                        lastBorderSegment = this.equalWithEpsilon(vz.y,yb);
                        vb = this.createVertex(lastBorderSegment ? vz.x : xr, yb);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        va = vb;
                        // fall through

                        // walk upward along right side
                        lastBorderSegment = this.equalWithEpsilon(vz.x,xr);
                        vb = this.createVertex(xr, lastBorderSegment ? vz.y : yt);
                        edge = this.createBorderEdge(cell.site, va, vb);
                        iLeft++;
                        halfedges.splice(iLeft, 0, this.createHalfedge(edge, cell.site, null));
                        nHalfedges++;
                        if ( lastBorderSegment ) { break; }
                        // fall through

                    default:
                        throw "Voronoi.closeCells() > this makes no sense!";
                    }
                }
            iLeft++;
            }
        cell.closeMe = false;
        }
    };

// ---------------------------------------------------------------------------
// Debugging helper
/*
Voronoi.prototype.dumpBeachline = function(y) {
    console.log('Voronoi.dumpBeachline(%f) > Beachsections, from left to right:', y);
    if ( !this.beachline ) {
        console.log('  None');
        }
    else {
        var bs = this.beachline.getFirst(this.beachline.root);
        while ( bs ) {
            console.log('  site %d: xl: %f, xr: %f', bs.site.voronoiId, this.leftBreakPoint(bs, y), this.rightBreakPoint(bs, y));
            bs = bs.rbNext;
            }
        }
    };
*/

// ---------------------------------------------------------------------------
// Helper: Quantize sites

// rhill 2013-10-12:
// This is to solve https://github.com/gorhill/Javascript-Voronoi/issues/15
// Since not all users will end up using the kind of coord values which would
// cause the issue to arise, I chose to let the user decide whether or not
// he should sanitize his coord values through this helper. This way, for
// those users who uses coord values which are known to be fine, no overhead is
// added.

Voronoi.prototype.quantizeSites = function(sites) {
    var ε = this.ε,
        n = sites.length,
        site;
    while ( n-- ) {
        site = sites[n];
        site.x = Math.floor(site.x / ε) * ε;
        site.y = Math.floor(site.y / ε) * ε;
        }
    };

// ---------------------------------------------------------------------------
// Helper: Recycle diagram: all vertex, edge and cell objects are
// "surrendered" to the Voronoi object for reuse.
// TODO: rhill-voronoi-core v2: more performance to be gained
// when I change the semantic of what is returned.

Voronoi.prototype.recycle = function(diagram) {
    if ( diagram ) {
        if ( diagram instanceof this.Diagram ) {
            this.toRecycle = diagram;
            }
        else {
            throw 'Voronoi.recycleDiagram() > Need a Diagram object.';
            }
        }
    };

// ---------------------------------------------------------------------------
// Top-level Fortune loop

// rhill 2011-05-19:
//   Voronoi sites are kept client-side now, to allow
//   user to freely modify content. At compute time,
//   *references* to sites are copied locally.

Voronoi.prototype.compute = function(sites, bbox) {
    // to measure execution time
    var startTime = new Date();

    // init internal state
    this.reset();

    // any diagram data available for recycling?
    // I do that here so that this is included in execution time
    if ( this.toRecycle ) {
        this.vertexJunkyard = this.vertexJunkyard.concat(this.toRecycle.vertices);
        this.edgeJunkyard = this.edgeJunkyard.concat(this.toRecycle.edges);
        this.cellJunkyard = this.cellJunkyard.concat(this.toRecycle.cells);
        this.toRecycle = null;
        }

    // Initialize site event queue
    var siteEvents = sites.slice(0);
    siteEvents.sort(function(a,b){
        var r = b.y - a.y;
        if (r) {return r;}
        return b.x - a.x;
        });

    // process queue
    var site = siteEvents.pop(),
        siteid = 0,
        xsitex, // to avoid duplicate sites
        xsitey,
        cells = this.cells,
        circle;

    // main loop
    for (;;) {
        // we need to figure whether we handle a site or circle event
        // for this we find out if there is a site event and it is
        // 'earlier' than the circle event
        circle = this.firstCircleEvent;

        // add beach section
        if (site && (!circle || site.y < circle.y || (site.y === circle.y && site.x < circle.x))) {
            // only if site is not a duplicate
            if (site.x !== xsitex || site.y !== xsitey) {
                // first create cell for new site
                cells[siteid] = this.createCell(site);
                site.voronoiId = siteid++;
                // then create a beachsection for that site
                this.addBeachsection(site);
                // remember last site coords to detect duplicate
                xsitey = site.y;
                xsitex = site.x;
                }
            site = siteEvents.pop();
            }

        // remove beach section
        else if (circle) {
            this.removeBeachsection(circle.arc);
            }

        // all done, quit
        else {
            break;
            }
        }

    // wrapping-up:
    //   connect dangling edges to bounding box
    //   cut edges as per bounding box
    //   discard edges completely outside bounding box
    //   discard edges which are point-like
    this.clipEdges(bbox);

    //   add missing edges in order to close opened cells
    this.closeCells(bbox);

    // to measure execution time
    var stopTime = new Date();

    // prepare return values
    var diagram = new this.Diagram();
    diagram.cells = this.cells;
    diagram.edges = this.edges;
    diagram.vertices = this.vertices;
    diagram.execTime = stopTime.getTime()-startTime.getTime();

    // clean up
    this.reset();

    return diagram;
    };

module.exports = Voronoi;
});

var areaPolygon = function (points,signed) {
  var l = points.length;
  var det = 0;
  var isSigned = signed || false;

  points = points.map(normalize);
  if (points[0] != points[points.length -1])  
    points = points.concat(points[0]);

  for (var i = 0; i < l; i++)
    det += points[i].x * points[i + 1].y
      - points[i].y * points[i + 1].x;
  if (isSigned)
    return det / 2
  else
    return Math.abs(det) / 2
};

function normalize(point) {
  if (Array.isArray(point))
    return {
      x: point[0],
      y: point[1]
    }
  else
    return point
}

var pointInPolygon = function (point, vs) {
    // ray-casting algorithm based on
    // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    
    var x = point[0], y = point[1];
    
    var inside = false;
    for (var i = 0, j = vs.length - 1; i < vs.length; j = i++) {
        var xi = vs[i][0], yi = vs[i][1];
        var xj = vs[j][0], yj = vs[j][1];
        
        var intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    
    return inside;
};

var classCallCheck = function (instance, Constructor) {
  if (!(instance instanceof Constructor)) {
    throw new TypeError("Cannot call a class as a function");
  }
};

var createClass = function () {
  function defineProperties(target, props) {
    for (var i = 0; i < props.length; i++) {
      var descriptor = props[i];
      descriptor.enumerable = descriptor.enumerable || false;
      descriptor.configurable = true;
      if ("value" in descriptor) descriptor.writable = true;
      Object.defineProperty(target, descriptor.key, descriptor);
    }
  }

  return function (Constructor, protoProps, staticProps) {
    if (protoProps) defineProperties(Constructor.prototype, protoProps);
    if (staticProps) defineProperties(Constructor, staticProps);
    return Constructor;
  };
}();

/**
 * VectorTools is not instanciable and provide some static functions for
 * computing things about vectors.
 */
var VectorTools = function () {
  function VectorTools() {
    classCallCheck(this, VectorTools);
  }

  createClass(VectorTools, null, [{
    key: "crossProduct",


    /**
     * performs a cross product between v1 and v2.
     * @param  {Array} v1 - a vector [x, y, z]. To use with 2D vectors, just use [x, y, 0]
     * @param  {Array} v2 - a vector [x, y, z]. To use with 2D vectors, just use [x, y, 0]
     * @param  {Boolean} normalize - will force normalization of the output vector if true (default: false)
     * @return {Array} a vector [x, y, z], result of the cross product
     */
    value: function crossProduct(v1, v2) {
      var normalize = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : false;

      var normalVector = [v1[1] * v2[2] - v1[2] * v2[1], (v1[0] * v2[2] - v1[2] * v2[0]) * -1, v1[0] * v2[1] - v1[1] * v2[0]];

      if (normalize) normalVector = VectorTools.normalize(normalVector);

      return normalVector;
    }

    /**
     * Perform the dot product of two vectors. They need to be of the same dimension
     * but they can be 2D, 3D or aother.
     * Note: If v1 and v2 are normalize, the dot product is also the cosine
     * @param  {Array} v1  - a vector [x, y] or [x, y, z]
     * @param  {Array} v2 - a vector [x, y] or [x, y, z]
     * @return {Number} the dot product
     */

  }, {
    key: "dotProduct",
    value: function dotProduct(v1, v2) {
      if (v1.length != v2.length) {
        console.log("ERROR: v1 and v2 must be the same size to compute a dot product");
        return null;
      }

      var sum = 0;

      for (var i = 0; i < v1.length; i++) {
        sum += v1[i] * v2[i];
      }

      return sum;
    }

    /**
     * Normalizes a 3D vector
     * @param  {Array} v - 3D vector to normalize
     * @return {Array} the normalized 3D vector
     */

  }, {
    key: "normalize",
    value: function normalize(v) {
      var n = VectorTools.getNorm(v);
      var normalizedV = [v[0] / n, v[1] / n, v[2] / n];
      return normalizedV;
    }

    /**
     * return the norm (length) of a vector [x, y, z]
     * @param  {Array} v - 3D vector to get the norm of
     * @return {Number} the norm
     */

  }, {
    key: "getNorm",
    value: function getNorm(v) {
      return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    /*
       Args:
        v: array - [x, y, z]
        m: array[array] -
            [[a, b, c],
             [d, e, f],
             [g, h, i]]
       Return rotated vector:
        [ax + by + cz,
         dx + ey + fz,
         gx + hy + iz]
    */

    /**
     * rotate the vector v using the rotation matrix m.
     * @param  {Array} v - 3D vector
     * @param  {Array} m  - matrix [[a, b, c], [d, e, f], [g, h, i]]
     * @return {Array} the rotated vector [ax + by + cz, dx + ey + fz, gx + hy + iz]
     */

  }, {
    key: "rotate",
    value: function rotate(v, m) {
      var vRot = [v[0] * m[0][0] + v[1] * m[0][1] + v[2] * m[0][2], v[0] * m[1][0] + v[1] * m[1][1] + v[2] * m[1][2], v[0] * m[2][0] + v[1] * m[2][1] + v[2] * m[2][2]];

      return vRot;
    }

    /**
     * Compute the angle p1p2p3 in radians.
     * Does not give the sign, just absolute angle!
     * @param  {Array} p1 - a 3D point [x, y, z]
     * @param  {Array} p2 - a 3D point [x, y, z]
     * @param  {Array} p3 - a 3D point [x, y, z]
     * @return {Number}    [description]
     */

  }, {
    key: "getAnglePoints",
    value: function getAnglePoints(p1, p2, p3) {

      //  the configuration is like that:
      //
      //  p1-----p2
      //        /
      //       /
      //      p3

      var v_p2p1 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];

      var v_p2p3 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];

      // normalizing those vectors
      v_p2p1 = VectorTools.normalize(v_p2p1);
      v_p2p3 = VectorTools.normalize(v_p2p3);

      var cosine = VectorTools.dotProduct(v_p2p1, v_p2p3);
      var angleRad = Math.acos(cosine);

      return angleRad;
    }

    /**
    * Checks if the 2D vector u crosses the 2D vector v. The vector u goes from point u1
    * to point u2 and vector v goes from point v1 to point v2.
    * Based on Gavin from SO http://bit.ly/2oNn741 reimplemented in JS
    * @param  {Array} u1 - first point of the u vector as [x, y]
    * @param  {Array} u2 - second point of the u vector as [x, y]
    * @param  {Array} v1 - first point of the v vector as [x, y]
    * @param  {Array} v2 - second point of the v vector as [x, y]
    * @return {Array|null} crossing point as [x, y] or null if vector don't cross.
    */

  }, {
    key: "vector2DCrossing",
    value: function vector2DCrossing(u1, u2, v1, v2) {
      var meet = null;
      var s1x = u2[0] - u1[0];
      var s1y = u2[1] - u1[1];
      var s2x = v2[0] - v1[0];
      var s2y = v2[1] - v1[1];

      var s = (-s1y * (u1[0] - v1[0]) + s1x * (u1[1] - v1[1])) / (-s2x * s1y + s1x * s2y);
      var t = (s2x * (u1[1] - v1[1]) - s2y * (u1[0] - v1[0])) / (-s2x * s1y + s1x * s2y);

      if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        meet = [u1[0] + t * s1x, u1[1] + t * s1y];
      }

      return meet;
    }

    /**
    * Get the distance between the 2D point p and a 2D segment u (defined by its points u1 and u2)
    * Stolen from SO http://bit.ly/2oQG3yX
    * @param  {Array} p  - a point as [x, y]
    * @param  {Array} u1 - first point of the u vector as [x, y]
    * @param  {Array} u2 - second point of the u vector as [x, y]
    * @return {Number} the distance
    */

  }, {
    key: "pointToSegmentDistance",
    value: function pointToSegmentDistance(p, u1, u2) {
      var A = p[0] - u1[0];
      var B = p[1] - u1[1];
      var C = u2[0] - u1[0];
      var D = u2[1] - u1[1];

      var dot = A * C + B * D;
      var lenSq = C * C + D * D;
      var param = -1;
      if (lenSq != 0) //in case of 0 length line
        param = dot / lenSq;

      var xx, yy;

      if (param < 0) {
        xx = u1[0];
        yy = u1[1];
      } else if (param > 1) {
        xx = u2[0];
        yy = u2[1];
      } else {
        xx = u1[0] + param * C;
        yy = u1[1] + param * D;
      }

      var dx = p[0] - xx;
      var dy = p[1] - yy;
      return Math.sqrt(dx * dx + dy * dy);
    }
  }]);
  return VectorTools;
}();

/**
 * ConvexPolygon is a simple approach of a polygon, it's actually a simple approcah
 * of what is a convex polygon and is mainly made to be used in the context of
 * a polygon representing a cell of a Voronoi diagram.
 * Here, a polygon is a list of 2D points that represent a convex polygon, with the
 * list of point being no closed (= the last point is not a repetition of the first).
 * The list of points describing a polygon can not be modified later.
 */

var ConvexPolygon = function () {

  /**
   * Contructor.
   * @param {Array} points - can be an array of [x, y] (both being of type Number)
   * or it can be an array of {x: Number, y: Number}. Depending on what is the source
   * of points, both exist.
   * Note: there is no integrity checking that the given list of point represent
   * an actually convex polygon.
   */
  function ConvexPolygon(points) {
    classCallCheck(this, ConvexPolygon);

    this._isValid = false;
    this._hull = null;
    this._area = null;

    if (Array.isArray(points) && points.length > 2) {
      var polygonPoints = null;
      // we are dealing with a [ {x, y}, {x, y}, ... ] polygon
      if ("x" in points[0] && "y" in points[0]) {
        polygonPoints = points.map(function (p) {
          return [p.x, p.y];
        });

        // we are dealing with a [ [x, y], [x, y], ... ]
      } else if (Array.isArray(points[0])) {
        polygonPoints = points.map(function (p) {
          return [p[0], p[1]];
        });
      }

      if (polygonPoints) {
        polygonPoints = ConvexPolygon.removeDuplicateVertices(polygonPoints);
        this._hull = ConvexPolygon.orderPolygonPoints(polygonPoints);
      }

      this._isValid = !!this._hull;
    }
  }

  /**
   * Get if the polygon is valid
   * @return {Boolean} [description]
   */


  createClass(ConvexPolygon, [{
    key: 'isValid',
    value: function isValid() {
      return this._isValid;
    }

    /**
     * [STATIC]
     * Removes the duplicates of a hull so that a polygon hull does not have twice
     * the same [x, y] points;
     * @param {Array} polygonPoints - an array of [x, y]
     * @return {Array} an array of [x, y], with duplicates removed
     */

  }, {
    key: 'getArea',


    /**
     * Get the area of a this polygon. If never computed, computes it and stores it
     * @return {Number} the area
     */
    value: function getArea() {
      if (!this._area) {
        this._area = areaPolygon(this._hull);
      }
      return this._area;
    }

    /**
    * This is reliable ONLY in the context of Voronoi cells! This is NOT a generally
    * valid way to find the intersection polygon between 2 convex polygons!
    * For this context only, we believe it's much faster tho.
    * @param {Polygon} anotherPolygon - another polygon the get the intersection with.
    */

  }, {
    key: 'getCellReplacementIntersection',
    value: function getCellReplacementIntersection(anotherPolygon) {
      var h1 = this.getHull();
      var h2 = anotherPolygon.getHull();
      var eps = 0.0001;

      var intersecVectices = [];
      var nbVerticesP1 = h1.length;
      var nbVerticesP2 = h2.length;

      // for each vertex of h1 ...
      for (var i = 0; i < nbVerticesP1; i++) {
        var u1 = h1[i];
        var u2 = h1[(i + 1) % nbVerticesP1];

        // case 1: all the vertices of h1 that are inside h2 have to be part of the list
        var inside = pointInPolygon(h1[i], h2);
        if (inside) intersecVectices.push(h1[i]);

        // for each vertex of h2 ...
        for (var j = 0; j < nbVerticesP2; j++) {
          // case 1 bis: all the vertices of h2 that are inside h1 have to be part of the list.
          // no need to run that every i loop
          if (i === 0) {
            var inside = pointInPolygon(h2[j], h1);
            if (inside) intersecVectices.push(h2[j]);
          }

          var v1 = h2[j];
          var v2 = h2[(j + 1) % nbVerticesP2];

          // case 2: get the intersection points between the edges of this poly and
          // the edges of anotherPolygon.
          var intersectPoint = VectorTools.vector2DCrossing(u1, u2, v1, v2);
          if (intersectPoint) {
            intersecVectices.push(intersectPoint);
          }

          // case 3: a vertex of a polygon is ON/ALONG an edge of the other polygon
          // note: this case can seem like "ho, that's an unfortunate minor case" but in the context
          // of voronoi cell replacement, this happens A LOT!

          // distance between the point v1 (that belongs to h2) and the edge u1u2 (that belongs to h1)
          var dv1u = VectorTools.pointToSegmentDistance(v1, u1, u2);
          if (dv1u < eps) {
            intersecVectices.push(v1);
          }

          // distance between the point u1 (that belongs to h1) and the edge v1v2 (that belongs to h2)
          var du1v = VectorTools.pointToSegmentDistance(u1, v1, v2);
          if (du1v < eps) {
            intersecVectices.push(u1);
          }
        }
      }

      var interPolygon = new ConvexPolygon(intersecVectices);
      return interPolygon;
    }

    /**
     * Get the Array of vertices. This is a reference and the point should not be
     * modified!
     * @return {Array} the points of the hull
     */

  }, {
    key: 'getHull',
    value: function getHull() {
      return this._hull;
    }

    /**
    * Compare with another polygon and tells if it's the same.
    * Predicate: polygons are convex + the first point of the list starts at noon and
    * the following are going clock-wise.
    * @param {ConvexPolygon} otherPolygon - another polygon
    * @return {Boolean} true is the same, false if not
    */

  }, {
    key: 'isSame',
    value: function isSame(otherPolygon) {
      var eps = 0.0001;
      var otherHull = otherPolygon.getHull();

      if (this._hull.length !== otherHull.length) return false;

      for (var i = 0; i < otherHull.length; i++) {
        if (Math.abs(otherHull[i][0] - this._hull[i][0]) > eps || Math.abs(otherHull[i][1] - this._hull[i][1]) > eps) return false;
      }

      return true;
    }
  }], [{
    key: 'removeDuplicateVertices',
    value: function removeDuplicateVertices(polygonPoints) {
      var newPolyPoints = [polygonPoints[0]];
      var eps = 0.0001;

      for (var i = 1; i < polygonPoints.length; i++) {
        var alreadyIn = false;
        for (var j = 0; j < newPolyPoints.length; j++) {
          var xDiff = Math.abs(polygonPoints[i][0] - newPolyPoints[j][0]);
          var yDiff = Math.abs(polygonPoints[i][1] - newPolyPoints[j][1]);
          //var zDiff = Math.abs(polygonPoints[i][2] - newPolyPoints[j][2]);

          if (xDiff < eps && yDiff < eps /*&& (zDiff < eps)*/) {
              alreadyIn = true;
              break;
            }
        }
        if (!alreadyIn) {
          newPolyPoints.push(polygonPoints[i]);
        }
      }
      return newPolyPoints;
    }

    /**
     * Get the center coordinate of a polygon by averaging all the pointd of the hull
     * @param {Array} polygonPoints - an array of [x, y]
     * @return {Array} the center as [x, y]
     */

  }, {
    key: 'getPolygonCenter',
    value: function getPolygonCenter(polygonPoints) {
      var nbVertice = polygonPoints.length;

      // find the center of the polygon
      var xAvg = 0;
      var yAvg = 0;
      var zAvg = 0;

      for (var v = 0; v < nbVertice; v++) {
        xAvg += polygonPoints[v][0];
        yAvg += polygonPoints[v][1];
        zAvg += polygonPoints[v][2];
      }

      xAvg /= nbVertice;
      yAvg /= nbVertice;
      zAvg /= nbVertice;
      var center = [xAvg, yAvg, zAvg];

      return center;
    }

    /**
     * A list of polygon points representing a convex hull are not always listed in
     * a clock-wise of ccw order, by default, we cannot count on it. This methods
     * reorder the vertices of the in a clock-wise fashion, starting by the one that
     * at noon or immediately after.
     * Note: this is necessary to compute the area and to compare two polygons.
     * @param {Array} polygonPoints - an array of [x, y]
     * @return {Array} an array of [x, y]
     */

  }, {
    key: 'orderPolygonPoints',
    value: function orderPolygonPoints(polygonPoints) {
      var nbVertice = polygonPoints.length;
      var center = ConvexPolygon.getPolygonCenter(polygonPoints);

      // for each, we have .vertice (a [x, y, z] array) and .angle (rad angle to planePolygonWithAngles[0])
      var planePolygonWithAngles = new Array(nbVertice);
      var verticalRay = [0, 1, 0];

      for (var v = 0; v < nbVertice; v++) {
        var currentRay = [center[0] - polygonPoints[v][0], center[1] - polygonPoints[v][1], 0];

        var currentRayNormalized = VectorTools.normalize(currentRay);
        var cos = VectorTools.dotProduct(verticalRay, currentRayNormalized);
        var angle = Math.acos(cos);
        var currentPolygonNormal = VectorTools.crossProduct(verticalRay, currentRayNormalized);
        var planeNormal = [0, 0, 1];
        var angleSign = VectorTools.dotProduct(currentPolygonNormal, planeNormal) > 0 ? 1 : -1;
        angle *= angleSign;

        // having only positive angles is a trick for ordering the vertices of a polygon
        // always the same way: first vertex of the list is at noon or the one just after
        // noon in a clock-wise direction. Then, all the other vertices are follwing in CW.
        // Then, it's easy to figure if 2 polygons are the same.
        if (angle < 0) {
          angle = Math.PI + (Math.PI + angle);
        }

        planePolygonWithAngles[v] = { vertex: polygonPoints[v], angle: angle };
      }

      // sort vertices based on their angle to [0]
      planePolygonWithAngles.sort(function (a, b) {
        return a.angle - b.angle;
      });

      // make a array of vertex only (ordered)
      var orderedVertice = [];
      for (var v = 0; v < nbVertice; v++) {
        orderedVertice.push(planePolygonWithAngles[v].vertex);
      }

      return orderedVertice;
    }
  }]);
  return ConvexPolygon;
}();

/**
* A Cell instance defines a polygon from a Voronoi diagram and its seed.
*
*/

var Cell = function () {

  /**
  * Constructor, from a list of points and a seed.
  * There is no integrity checking to make sure the seed is actually within the polygon.
  * The polygon, since it's supposed to come from a Voronoi cell, is expected to be convex.
  * @param {Array} contourPoints - Array of points
  */
  function Cell(contourPoints, seed) {
    classCallCheck(this, Cell);

    this._hash = Cell.genarateHash(seed.x, seed.y);
    this._polygon = new ConvexPolygon(contourPoints);
    this._seed = seed;
    this._isValid = this._polygon.isValid() && !!this._hash;
  }

  /**
  * Return if this cell is valid
  * @return {Boolean} true if valid, false if not
  */


  createClass(Cell, [{
    key: "isValid",
    value: function isValid() {
      return this._isValid;
    }

    /**
    * Get the hash of this cell
    * @return {String}
    */

  }, {
    key: "getHash",
    value: function getHash() {
      return this._hash;
    }

    /**
    * Get the seed of the cell.
    * @return {Object} should be of form {x: Number, y: Number, seedIndex: Number}
    */

  }, {
    key: "getSeed",
    value: function getSeed() {
      return this._seed;
    }

    /**
    * Get the polygon of this cell
    * @return {ConvexPolygon}
    */

  }, {
    key: "getPolygon",
    value: function getPolygon() {
      return this._polygon;
    }

    /**
    * [STATIC]
    * Get a unique hash from 2 floats. I am pretty sure these super simple stringified
    * coordinate hashes are not the most robust, though the point is mainly to be as
    * fast as possible.
    * @param {Number} n1 - a number, most likely the x of a coord
    * @param {Number} n2 - a number, most likely the y of a coord
    * @param {String} a (probably unique enough) hash
    */

  }, {
    key: "getArea",


    /**
    * Get the area of the cell (calls the polygon method)
    * @return {Number} the area
    */
    value: function getArea() {
      return this._polygon.getArea();
    }

    /**
    * Get the polygon resulting from the intersection of a voronoi cell with another
    * cell (this second cell comes from another voronoi diagramm where a seed has been added).
    * NOTE: This does NOT implement a standard polygon intersection algorithm, its use
    * is stricly for the use case of voronoi cells being replaced, making it quite faster
    * but not suitable for other cases.
    * @param {Cell} anotherCell - another cell to intersect with
    * @return {ConvexPolygon} the polygon resulting from the intersection
    */

  }, {
    key: "intersectWithCell",
    value: function intersectWithCell(anotherCell) {
      return this._polygon.getCellReplacementIntersection(anotherCell.getPolygon());
    }

    /**
    * Compare if this cell has the same polygon as another cell.
    * Read the doc of ConvexPolygon.getPolygon() for more info.
    * @param {Cell} anotherCell - another cell to compare the polygon with.
    * @return {Boolean} true if the same polygon, false if not
    */

  }, {
    key: "hasSamePolygon",
    value: function hasSamePolygon(otherCell) {
      return this._polygon.isSame(otherCell.getPolygon());
    }
  }], [{
    key: "genarateHash",
    value: function genarateHash(n1, n2) {
      //var hash = (i1 < base ? "0" : '') + i1.toString(16) +
      //           (i2 < base ? "0" : '') + i2.toString(16)

      //
      var hash = n1.toString() + "_" + n2.toString();
      return hash;
    }
  }]);
  return Cell;
}();

/**
* A CellCollection stores Cell objects. It is the interface between the Voronoi diagram
* and the Cells/Polygons.
*/

var CellCollection = function () {

  /**
  * The constructor does not take any param.
  */
  function CellCollection() {
    classCallCheck(this, CellCollection);

    this._cells = {};
    this._cellArray = [];
    this._addedSeedHash = null;
  }

  /**
  * Builds the collection of cell using a voronoi diagram object from the Javascript-Voronoi
  * library by Gorhill ( http://bit.ly/2FoapPI )
  */


  createClass(CellCollection, [{
    key: 'buildFromVoronoiDiagram',
    value: function buildFromVoronoiDiagram(vd) {
      var vCells = vd.cells;

      // for each cell
      for (var i = 0; i < vCells.length; i++) {
        var cellhes = vCells[i].halfedges;
        var points = [];

        for (var j = 0; j < cellhes.length; j++) {
          points.push(cellhes[j].edge.va);
          points.push(cellhes[j].edge.vb);
        }

        // build a cell
        var cell = new Cell(points, vCells[i].site);
        if (cell.isValid()) {
          this._cells[cell.getHash()] = cell;
          this._cellArray.push(cell);
        }
      }
    }

    /**
    * Get the list of hashes (ID) of cells in the collection
    * @return {Array} all the hashes
    */

  }, {
    key: 'getCellHashes',
    value: function getCellHashes() {
      return Object.keys(this._cells);
    }

    /**
    * Get the cell of _this_ collection that has this hash
    * @param {String} hash - the unique hash of the cell
    * @return {Cell} a Cell instance of null if hash not found
    */

  }, {
    key: 'getCell',
    value: function getCell(hash) {
      if (hash in this._cells) return this._cells[hash];else return null;
    }

    /**
    * Store the hash of a given position of a seed, meaning, store a reference to a cell.
    * This seed is the one added to the 'another' diagram, when we introduce a pixel as a seed.
    * For instance, the seed-only cell collection does not have a reference to a seed.
    * @param {Number} x - the x position of the seed to reference
    * @param {Number} y - the y position of the seed to reference
    */

  }, {
    key: 'referenceSeed',
    value: function referenceSeed(x, y) {
      this._addedSeedHash = Cell.genarateHash(x, y);
    }

    /**
    * Get the cell that was referenced by referenceSeed()
    * @return {Cell}
    */

  }, {
    key: 'getReferencedSeedCell',
    value: function getReferencedSeedCell() {
      return this._cells[this._addedSeedHash];
    }

    /**
    * When comparing a seed-only CellCollection with a seed-and-one-pixel CellCollection
    * with getStolenAreaInfo(), we get a list of original cell index (from the seed-only collection)
    * as well as the area ratio that comes from them to build the pixel-based cell
    * @param {CellCollection} anotherCellCollection - a cell collection with the original seeds + 1 pixel as a seed
    * @return {Array} list of {seedIndex: Number, weight: Number}
    */

  }, {
    key: 'getStolenAreaInfo',
    value: function getStolenAreaInfo(anotherCellCollection) {
      var addedCell = anotherCellCollection.getReferencedSeedCell();
      var addedCellArea = addedCell.getArea();
      var modifiedCells = [];

      for (var hash in this._cells) {
        if (hash === addedCell.getHash()) continue;

        var originalCell = this._cells[hash];
        var newCell = anotherCellCollection.getCell(hash);
        var isSame = originalCell.hasSamePolygon(newCell);

        if (!isSame) {
          var stolenRatio = 0;
          var stolenPoly = addedCell.intersectWithCell(originalCell);

          if (stolenPoly.isValid()) {
            var stolenArea = stolenPoly.getArea();
            stolenRatio = stolenArea / addedCell.getArea();
          }

          stolenRatio = Math.round(stolenRatio * 10000) / 10000;

          //modifiedCells.push( {cell: originalCell, stolenRatio: stolenRatio} );
          modifiedCells.push({ seedIndex: originalCell.getSeed().seedIndex, weight: stolenRatio });
        }
      }

      return modifiedCells;
    }
  }]);
  return CellCollection;
}();

/**
* The Interpolator is the API provider of natninter.
*/

var Interpolator = function () {

  /**
  * No param constructor
  */
  function Interpolator() {
    classCallCheck(this, Interpolator);

    this._seeds = [];
    this._output = { width: 256, height: 256 };
    this._samplingMap = [];
    this._recomputeMap = true;

    // voronoi diagram of seeds only
    this._seedCellCollection = null;

    this._onProgressCallback = null;
  }

  /**
   * Define a callback for the progress of making the map and the progress of making the output image
   * It will be called with 2 args: the progress in percentage (Number), and a description string
   * @param  {Function} cb - callback function
   */


  createClass(Interpolator, [{
    key: 'onProgress',
    value: function onProgress(cb) {
      if (cb && typeof cb === "function") {
        this._onProgressCallback = cb;
      }
    }

    /**
    * Removes all the seeds
    */

  }, {
    key: 'cleanSeeds',
    value: function cleanSeeds() {
      this._seeds = [];
      this._recomputeMap = true;
    }

    /**
    * Add a seed to the interpolator. Note that the seed is not copied, only its reference is.
    * This is convenient for updating the value of the seed from the outside and ecompute
    * the interpolation without having to recompute the whole weight map.
    * @param {Object} seed - of form {x: Number, y: Number, value: Number}
    */

  }, {
    key: 'addSeed',
    value: function addSeed(seed) {
      if ("x" in seed && "y" in seed && "value" in seed) {
        this._seeds.push(seed);
        this._recomputeMap = true;
      }
    }

    /**
    * Add an array of seed
    * @param {Array} seedArr - array of seeds, where each seed is of form {x: Number, y: Number, value: Number}
    */

  }, {
    key: 'addSeeds',
    value: function addSeeds(seedArr) {
      for (var i = 0; i < seedArr.length; i++) {
        this.addSeed(seedArr[i]);
      }
    }

    /**
     * Get the number of seeds, mainly to add a validation steop after adding them.
     * @return {Number} The number of seed
     */

  }, {
    key: 'getNumberOfSeeds',
    value: function getNumberOfSeeds() {
      return this._seeds.length;
    }

    /**
    * Define the size of the output image
    * @param {Number} width - width of the output image
    * @param {Number} height - height of the output image
    */

  }, {
    key: 'setOutputSize',
    value: function setOutputSize(width, height) {
      if (width > 0 && height > 0) {
        this._output.width = width;
        this._output.height = height;
      }
    }

    /**
    * Check if all the seeds are inside the output area defined with `setOutputSize()`
    * @return {Boolean} true if all are inside, false if at leas one is out
    */

  }, {
    key: 'hasAllSeedsInside',
    value: function hasAllSeedsInside() {
      var size = this._output;
      return this._seeds.every(function (s) {
        return s.x >= 0 && s.y >= 0 && s.x < size.width && s.y < size.width;
      });
    }

    /**
    * Compute the sampling map. Automatically called by the update method when the
    * map was never computed or when a seed have been added since.
    * Though this method is not private and can be called to force recomputing the map.
    * Note: if you have already computed a map for the exact same seed positions and output size,
    * you can avoid recomputing it and use `getMap` and `setMap`.
    * @return {Boolean} true if the process went well, false if error generating the map
    */

  }, {
    key: 'generateMap',
    value: function generateMap() {
      if (!this.hasAllSeedsInside()) {
        console.log('ERR: some seeds are outside of the image. Use .setOutputSize() to change the output size or modify the seed.');
        return false;
      }
      this._generateSeedCells();

      var w = this._output.width;
      var h = this._output.height;
      this._samplingMap = new Array(w * h);

      // for each pixel of the output image
      for (var i = 0; i < w; i++) {
        if (this._onProgressCallback) this._onProgressCallback(Math.round((i + 1) / w * 100), "sampling map");

        for (var j = 0; j < h; j++) {

          var seedIndex = this.isAtSeedPosition(i, j);
          var index1D = i + j * w;

          if (seedIndex === -1) {
            //this._samplingMap[ index1D ] = this._generatePixelCells(i, j);
            var pixelCellCollection = this._generatePixelCells(i, j);
            var stolenAreaData = this._seedCellCollection.getStolenAreaInfo(pixelCellCollection);
            this._samplingMap[index1D] = stolenAreaData;
          } else {
            this._samplingMap[index1D] = [{ seedIndex: seedIndex, weight: 1 }];
          }
        }
      }

      return true;
    }

    /**
    * Get the sampling map object
    */

  }, {
    key: 'getMap',
    value: function getMap(m) {
      return this._samplingMap;
    }

    /**
    * When you don't want to recompute the sampling map with `computeMap()` and reuse
    * the exact same seed position and output size (of course, seed values can change)
    * @param {Array} map - an already existing sampling map
    */

  }, {
    key: 'setMap',
    value: function setMap(map) {
      if (map.length === this._output.width * this._output.height) {
        this._samplingMap = map;
        return true;
      } else {
        console.log("The sampling map must be an 1D array of size output.x*output.y");
        return false;
      }
    }

    /**
    * is the given position at the position of a seed?
    * @param {Number} i - position along x axis
    * @param {Number} j - position along y axis
    * @return {Boolean} -1 if not, of the index of the seed if yes
    */

  }, {
    key: 'isAtSeedPosition',
    value: function isAtSeedPosition(i, j) {
      return this._seeds.findIndex(function (elem) {
        return elem.x == i && elem.y == j;
      });
    }

    /**
    * [PRIVATE]
    * Generate the voronoi diagram where sites are only seeds
    */

  }, {
    key: '_generateSeedCells',
    value: function _generateSeedCells() {
      var bbox = { xl: 0, xr: this._output.width, yt: 0, yb: this._output.height };
      var sites = this._seeds.map(function (s, i) {
        return { x: s.x, y: s.y, seedIndex: i };
      });

      var voronoi = new rhillVoronoiCore();
      var seedVoronoiDiagram = voronoi.compute(sites, bbox);

      this._seedCellCollection = new CellCollection();
      this._seedCellCollection.buildFromVoronoiDiagram(seedVoronoiDiagram);
    }
  }, {
    key: '_generatePixelCells',
    value: function _generatePixelCells(i, j) {
      var voronoi = new rhillVoronoiCore();
      var bbox = { xl: 0, xr: this._output.width, yt: 0, yb: this._output.height };
      var sites = this._seeds.map(function (s, i) {
        return { x: s.x, y: s.y, seedIndex: i };
      });
      sites.push({ x: i, y: j, seedIndex: -1 });
      var pixelVoronoiDiagram = voronoi.compute(sites, bbox);

      var pixelCellCollection = new CellCollection();
      pixelCellCollection.buildFromVoronoiDiagram(pixelVoronoiDiagram);
      pixelCellCollection.referenceSeed(i, j);

      return pixelCellCollection;
    }

    /**
    * Generate the output image as a floating points 1D array representing a 2D (1band)
    * image.
    * @return {Object} of form {_data: Float32Array, _metadata: {width: Number, height: Number}}
    */

  }, {
    key: 'generateImage',
    value: function generateImage() {
      var l = this._output.width * this._output.height;
      var outImg = new Float32Array(l);
      var seeds = this._seeds;
      var map = this._samplingMap;

      for (var i = 0; i < l; i++) {
        if (this._onProgressCallback) this._onProgressCallback(Math.round((i + 1) / l * 100), "output image");

        var pixelMap = map[i];
        var sum = 0;

        for (var m = 0; m < pixelMap.length; m++) {
          sum += pixelMap[m].weight * seeds[pixelMap[m].seedIndex].value;
        }

        outImg[i] = sum;
      }

      return {
        _metadata: {
          width: this._output.width,
          height: this._output.height
        },
        _data: outImg
      };
    }
  }]);
  return Interpolator;
}();

exports.Interpolator = Interpolator;

Object.defineProperty(exports, '__esModule', { value: true });

})));
//# sourceMappingURL=natninter.umd.js.map
