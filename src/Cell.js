import { ConvexPolygon } from "./ConvexPolygon.js";


/**
* A Cell instance defines a polygon from a Voronoi diagram and its seed.
*/
class Cell {

  /**
  * Build a cell.
  * @param {Array}
  */
  constructor( contourPoints, seed ){
    this._hash = Cell.getHash( seed.x, seed.y )
    this._polygon= new ConvexPolygon( contourPoints );
    this._seed = seed;
    this._isValid = (this._polygon.isValid() && !!this._hash);
  }


  isValid(){
    return this._isValid;
  }


  getHash(){
    return this._hash;
  }


  getSeed(){
    return this._seed;
  }


  getPolygon(){
    return this._polygon;
  }


  // get a unique hash from 2 floats
  static getHash( i1, i2 ){
    //return md5( new Float32Array([f1, f2]) );
    var base = 16;
    var hash = (i1 < base ? "0" : '') + i1.toString(16) +
               (i2 < base ? "0" : '') + i2.toString(16)
    return hash;
  }


  getArea(){
    return this._polygon.getArea();
  }


  intersectWithCell( anotherCell ){
    return this._polygon.getIntersection( anotherCell.getPolygon() );
  }


  hasSamePolygon( otherCell ){
    return this._polygon.isSame( otherCell.getPolygon() );
  }



}

export { Cell };