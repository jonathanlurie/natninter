import AreaPolygon from "area-polygon";
import greinerHormann from 'greiner-hormann';
import { VectorTools } from './VectorTools.js'


class ConvexPolygon {

  constructor( points ){
    this._isValid = false;
    this._hull = null;

    if( Array.isArray( points ) && points.length > 2 ){
      var polygonPoints = null;
      // we are dealing with a [ {x, y}, {x, y}, ... ] polygon
      if( "x" in points[0] && "y" in points[0] ){
        polygonPoints = points.map( function(p){ return [p.x, p.y, 0]} );

      // we are dealing with a [ [x, y], [x, y], ... ]
      }else if( Array.isArray( points[0] ) ){
        polygonPoints = points.map( function(p){ return [p[0], p[1], 0]} );
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
        var zDiff = Math.abs(polygonPoints[i][2] - newPolyPoints[j][2]);

        if( (xDiff < eps) && (yDiff < eps) && (zDiff < eps)){
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
    var nbVertice = polygonPoints.length;
    var center = ConvexPolygon.getPolygonCenter( polygonPoints );
    //console.log( center );

    // create normalized vectors from center to each vertex of the polygon
    var normalizedRays = [];
    for(var v=0; v<nbVertice; v++){
      var currentRay = [
        center[0] - polygonPoints[v][0],
        center[1] - polygonPoints[v][1],
        center[2] - polygonPoints[v][2]
      ];

      normalizedRays.push(VectorTools.normalize(currentRay));
    }

    // for each, we have .vertice (a [x, y, z] array) and .angle (rad angle to planePolygonWithAngles[0])
    var planePolygonWithAngles = [];

    // find the angle of each towards the first vertex
    planePolygonWithAngles.push({vertex: polygonPoints[0], angle: 0})
    for(var v=1; v<nbVertice; v++){
      var cos = VectorTools.dotProduct(normalizedRays[0], normalizedRays[v]);
      var angle = Math.acos(cos);
      var currentPolygonNormal = VectorTools.crossProduct(normalizedRays[0], normalizedRays[v], false);
      var planeNormal = [0, 0, 1];
      var angleSign = VectorTools.dotProduct(currentPolygonNormal, planeNormal)>0? 1:-1;
      angle *= angleSign;

      planePolygonWithAngles.push({vertex: polygonPoints[v], angle: angle})
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
    return orderedVertice;
  }



  getArea(){
    return AreaPolygon( this._hull );
  }


  getIntersection( anotherPolygon ){
    // TODO the intersection gives the input
    var vertices = greinerHormann.intersection( this.get2DHull() , anotherPolygon.get2DHull() );
    var interPolygon = new ConvexPolygon( vertices[0] );
    return interPolygon;
  }


  getHull(){
    return this._hull;
  }


  get2DHull(){
    return this._hull.map(function(p){ return [p[0], p[1]]});
  }

  isSame( otherPolygon ){
    var otherHull = otherPolygon.getHull();

    if( this._hull.length !== otherHull.length )
      return false;

    var isSame = true;

    for(var i=0; i<otherHull.length; i++){
      isSame = isSame || (otherHull[0]===this._hull[0] && otherHull[1]===this._hull[1])
      if( !isSame )
        break;
    }

    return isSame;
  }







}

export { ConvexPolygon };
