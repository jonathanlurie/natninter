import AreaPolygon from "area-polygon";
import PointInPolygon from 'point-in-polygon';

import { VectorTools } from './VectorTools.js'

/**
 * ConvexPolygon is a simple approach of a polygon, it's actually a simple approcah
 * of what is a convex polygon and is mainly made to be used in the context of
 * a polygon representing a cell of a Voronoi diagram.
 * Here, a polygon is a list of 2D points that represent a convex polygon, with the
 * list of point being no closed (= the last point is not a repetition of the first).
 * The list of points describing a polygon can not be modified later.
 */
class ConvexPolygon {

  /**
   * Contructor.
   * @param {Array} points - can be an array of [x, y] (both being of type Number)
   * or it can be an array of {x: Number, y: Number}. Depending on what is the source
   * of points, both exist.
   * Note: there is no integrity checking that the given list of point represent
   * an actually convex polygon.
   */
  constructor( points ){
    this._isValid = false;
    this._hull = null;
    this._area = null;

    if( Array.isArray( points ) && points.length > 2 ){
      var polygonPoints = null;
      // we are dealing with a [ {x, y}, {x, y}, ... ] polygon
      if( "x" in points[0] && "y" in points[0] ){
        polygonPoints = points.map( function(p){ return [p.x, p.y]} );

      // we are dealing with a [ [x, y], [x, y], ... ]
      }else if( Array.isArray( points[0] ) ){
        polygonPoints = points.map( function(p){ return [p[0], p[1]]} );
      }

      if( polygonPoints ){
        polygonPoints = ConvexPolygon.removeDuplicateVertices( polygonPoints )
        this._hull = ConvexPolygon.orderPolygonPoints( polygonPoints );
      }

      this._isValid = !!this._hull
    }
  }


  /**
   * Get if the polygon is valid
   * @return {Boolean} [description]
   */
  isValid(){
    return this._isValid;
  }


  /**
   * [STATIC]
   * Removes the duplicates of a hull so that a polygon hull does not have twice
   * the same [x, y] points;
   * @param {Array} polygonPoints - an array of [x, y]
   * @return {Array} an array of [x, y], with duplicates removed
   */
  static removeDuplicateVertices( polygonPoints ){
    var newPolyPoints = [ polygonPoints[0] ];
    var eps = 0.0001;

    for(var i=1; i<polygonPoints.length; i++){
      let alreadyIn = false;
      for(var j=0; j<newPolyPoints.length; j++){
        var xDiff = Math.abs(polygonPoints[i][0] - newPolyPoints[j][0]);
        var yDiff = Math.abs(polygonPoints[i][1] - newPolyPoints[j][1]);
        //var zDiff = Math.abs(polygonPoints[i][2] - newPolyPoints[j][2]);

        if( (xDiff < eps) && (yDiff < eps) /*&& (zDiff < eps)*/){
          alreadyIn = true;
          break
        }
      }
      if( !alreadyIn ){
        newPolyPoints.push( polygonPoints[i] );
      }
    }
    return newPolyPoints;
  }


  /**
   * Get the center coordinate of a polygon by averaging all the pointd of the hull
   * @param {Array} polygonPoints - an array of [x, y]
   * @return {Array} the center as [x, y]
   */
  static getPolygonCenter( polygonPoints ){
    var nbVertice = polygonPoints.length;

    // find the center of the polygon
    var xAvg = 0;
    var yAvg = 0;
    var zAvg = 0;

    for(var v=0; v<nbVertice; v++){
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
  static orderPolygonPoints( polygonPoints ){
    var nbVertice = polygonPoints.length;
    var center = ConvexPolygon.getPolygonCenter( polygonPoints );

    // for each, we have .vertice (a [x, y, z] array) and .angle (rad angle to planePolygonWithAngles[0])
    var planePolygonWithAngles = new Array( nbVertice );
    var verticalRay = [0, 1, 0];

    for(var v=0; v<nbVertice; v++){
      var currentRay = [
        center[0] - polygonPoints[v][0],
        center[1] - polygonPoints[v][1],
        0
      ];

      var currentRayNormalized = VectorTools.normalize(currentRay);
      var cos = VectorTools.dotProduct(verticalRay, currentRayNormalized);
      var angle = Math.acos(cos);
      var currentPolygonNormal = VectorTools.crossProduct(verticalRay, currentRayNormalized);
      var planeNormal = [0, 0, 1];
      var angleSign = VectorTools.dotProduct(currentPolygonNormal, planeNormal)>0? 1:-1;
      angle *= angleSign;

      // having only positive angles is a trick for ordering the vertices of a polygon
      // always the same way: first vertex of the list is at noon or the one just after
      // noon in a clock-wise direction. Then, all the other vertices are follwing in CW.
      // Then, it's easy to figure if 2 polygons are the same.
      if( angle < 0 ){
        angle = Math.PI + ( Math.PI + angle );
      }

      planePolygonWithAngles[v] = {vertex: polygonPoints[v], angle: angle};
    }

    // sort vertices based on their angle to [0]
    planePolygonWithAngles.sort(function(a, b){
      return a.angle - b.angle;
    });

    // make a array of vertex only (ordered)
    var orderedVertice = [];
    for(var v=0; v<nbVertice; v++){
      orderedVertice.push(planePolygonWithAngles[v].vertex);
    }

    return orderedVertice;
  }



  /**
   * Get the area of a this polygon. If never computed, computes it and stores it
   * @return {Number} the area
   */
  getArea(){
    if (! this._area){
      this._area = AreaPolygon( this._hull );
    }
    return this._area;
  }


  /**
  * This is reliable ONLY in the context of Voronoi cells! This is NOT a generally
  * valid way to find the intersection polygon between 2 convex polygons!
  * For this context only, we believe it's much faster tho.
  * @param {Polygon} anotherPolygon - another polygon the get the intersection with.
  */
  getCellReplacementIntersection( anotherPolygon ){
    var h1 = this.getHull();
    var h2 = anotherPolygon.getHull();
    var eps = 0.0001;


    var intersecVectices = [];
    var nbVerticesP1 = h1.length;
    var nbVerticesP2 = h2.length;

    // for each vertex of h1 ...
    for(var i=0; i<nbVerticesP1; i++){
      var u1 = h1[ i ];
      var u2 = h1[ (i+1)%nbVerticesP1 ];

      // case 1: all the vertices of h1 that are inside h2 have to be part of the list
      var inside = PointInPolygon( h1[ i ], h2);
      if( inside )
        intersecVectices.push( h1[ i ] );

      // for each vertex of h2 ...
      for(var j=0; j<nbVerticesP2; j++){
        // case 1 bis: all the vertices of h2 that are inside h1 have to be part of the list.
        // no need to run that every i loop
        if( i === 0 ){
          var inside = PointInPolygon( h2[ j ], h1);
          if( inside )
            intersecVectices.push( h2[ j ] );
        }

        var v1 = h2[ j ];
        var v2 = h2[ (j+1)%nbVerticesP2 ];

        // case 2: get the intersection points between the edges of this poly and
        // the edges of anotherPolygon.
        var intersectPoint = VectorTools.vector2DCrossing(u1, u2, v1, v2);
        if( intersectPoint ){
          intersecVectices.push( intersectPoint );
        }

        // case 3: a vertex of a polygon is ON/ALONG an edge of the other polygon
        // note: this case can seem like "ho, that's an unfortunate minor case" but in the context
        // of voronoi cell replacement, this happens A LOT!

        // distance between the point v1 (that belongs to h2) and the edge u1u2 (that belongs to h1)
        var dv1u = VectorTools.pointToSegmentDistance(v1, u1, u2);
        if( dv1u < eps ){
          intersecVectices.push( v1 );
        }

        // distance between the point u1 (that belongs to h1) and the edge v1v2 (that belongs to h2)
        var du1v = VectorTools.pointToSegmentDistance(u1, v1, v2);
        if( du1v < eps ){
          intersecVectices.push( u1 );
        }


      }
    }

    var interPolygon = new ConvexPolygon( intersecVectices );
    return interPolygon;
  }


  /**
   * Get the Array of vertices. This is a reference and the point should not be
   * modified!
   * @return {Array} the points of the hull
   */
  getHull(){
    return this._hull;
  }


  /**
  * Compare with another polygon and tells if it's the same.
  * Predicate: polygons are convex + the first point of the list starts at noon and
  * the following are going clock-wise.
  * @param {ConvexPolygon} otherPolygon - another polygon
  * @return {Boolean} true is the same, false if not
  */
  isSame( otherPolygon ){
    var eps = 0.0001;
    var otherHull = otherPolygon.getHull();

    if( this._hull.length !== otherHull.length )
      return false;

    for(var i=0; i<otherHull.length; i++){
      if( (Math.abs(otherHull[i][0] - this._hull[i][0])>eps || Math.abs(otherHull[i][1] - this._hull[i][1]) > eps ))
        return false;
    }

    return true;
  }

}

export { ConvexPolygon };
