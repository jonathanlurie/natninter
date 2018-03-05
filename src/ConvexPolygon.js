import AreaPolygon from "area-polygon";
//import greinerHormann from 'greiner-hormann';
import PolyBool from 'polybooljs';
import { VectorTools } from './VectorTools.js'


class ConvexPolygon {

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

  isValid(){
    return this._isValid;
  }


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

        /*
        if( (polygonPoints[i][0] - newPolyPoints[j][0]) < eps &&
            (polygonPoints[i][1] - newPolyPoints[j][1]) < eps &&
            (polygonPoints[i][2] - newPolyPoints[j][2]) < eps )
        {
          alreadyIn = true;
          break
        }
        */
      }
      if( !alreadyIn ){
        newPolyPoints.push( polygonPoints[i] );
      }
    }

    return newPolyPoints;
  }


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





  static orderPolygonPoints( polygonPoints ){
    var polygonPoints3D = polygonPoints.map(function(p){return [p[0], p[1], 0]});
    var nbVertice = polygonPoints3D.length;
    var center = ConvexPolygon.getPolygonCenter( polygonPoints3D );
    //console.log( center );

    // create normalized vectors from center to each vertex of the polygon
    var normalizedRays = [];
    for(var v=0; v<nbVertice; v++){
      var currentRay = [
        center[0] - polygonPoints3D[v][0],
        center[1] - polygonPoints3D[v][1],
        center[2] - polygonPoints3D[v][2]
      ];

      normalizedRays.push(VectorTools.normalize(currentRay));
    }

    // for each, we have .vertice (a [x, y, z] array) and .angle (rad angle to planePolygonWithAngles[0])
    var planePolygonWithAngles = [];
    var verticalRay = [0, 1, 0];

    // find the angle of each towards the first vertex
    //planePolygonWithAngles.push({vertex: polygonPoints3D[0], angle: 0})

    //for(var v=1; v<nbVertice; v++){
    for(var v=0; v<nbVertice; v++){
      var cos = VectorTools.dotProduct(verticalRay, normalizedRays[v]);
      var angle = Math.acos(cos);
      var currentPolygonNormal = VectorTools.crossProduct(verticalRay, normalizedRays[v], false);
      var planeNormal = [0, 0, 1];
      var angleSign = VectorTools.dotProduct(currentPolygonNormal, planeNormal)>0? 1:-1;
      angle *= angleSign;

      if( angle < 0 ){
        angle = Math.PI + ( Math.PI + angle );
      }

      planePolygonWithAngles.push({vertex: polygonPoints3D[v], angle: angle})
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

    // attribute the ordered array to this.planePolygo
    //newPolyPoints = newPolyPoints.map(function(p){return [p[0], p[1]]});
    orderedVertice = orderedVertice.map(function(p){return [p[0], p[1]]});
    return orderedVertice;
  }



  getArea(){
    if (! this._area){
      this._area = AreaPolygon( this._hull );
    }
    return this._area;
  }


  getIntersection( anotherPolygon ){
    // TODO the intersection gives the input
    var polyA = {
                regions: [
                  this.get2DHull()
                ],
                inverted: false
              }

    var polyB = {
                regions: [
                  anotherPolygon.get2DHull()
                ],
                inverted: false
              }

    var vertices = PolyBool.intersect( polyA , polyB );
    var interPolygon = new ConvexPolygon( vertices.regions[0] );
    return interPolygon;
  }


  getHull(){
    return this._hull;
  }


  get2DHull(){
    return this._hull.map(function(p){ return [p[0], p[1]]});
  }

  isSame( otherPolygon ){
    var eps = 0.0001;
    var otherHull = otherPolygon.getHull();

    if( this._hull.length !== otherHull.length )
      return false;

    for(var i=0; i<otherHull.length; i++){
      //if( (otherHull[i][0]!==this._hull[i][0] || otherHull[i][1]!==this._hull[i][1]))
      //  return false;

      if( (Math.abs(otherHull[i][0] - this._hull[i][0])>eps || Math.abs(otherHull[i][1] - this._hull[i][1]) > eps ))
        return false;
    }

    return true;
  }





}

export { ConvexPolygon };
