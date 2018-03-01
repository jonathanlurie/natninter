
import md5 from 'md5';
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
    this._polygon= new ConvexPolygon( contourPoints );
    this._hash = md5( new Float32Array([seed.x, seed.y]) );
    this._seed = seed;
    this._isValid = (this._polygon.isValid() && !!this._hash);
  }


  isValid(){
    return this._isValid;
  }

  getHash(){
    return this._hash;
  }



}

export { Cell };
