import convexHull from "convex-hull";

class ConvexPolygon {

  constructor( points ){
    this._isValid = false;
    this._hull = null;

    if( Array.isArray( points ) && points.length > 2 ){
      // we are dealing with a [ {x, y}, {x, y}, ... ] polygon
      if( "x" in points[0] && "y" in points[0] ){
        this._hull = convexHull( points.map( function(p){ return [p.x, p.y]} ) );

      // we are dealing with a [ [x, y], [x, y], ... ]
      }else if( Array.isArray( points[0] ) ){
        this._hull = convexHull( points )
      }

      this._isValid = !!this._hull
    }
  }

  isValid(){
    return this._isValid;
  }

}

export { ConvexPolygon };
