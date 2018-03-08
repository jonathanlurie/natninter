import { ConvexPolygon } from "./ConvexPolygon.js";


/**
* A Cell instance defines a polygon from a Voronoi diagram and its seed.
*
*/
class Cell {

  /**
  * Constructor, from a list of points and a seed.
  * There is no integrity checking to make sure the seed is actually within the polygon.
  * The polygon, since it's supposed to come from a Voronoi cell, is expected to be convex.
  * @param {Array} contourPoints - Array of points
  */
  constructor( contourPoints, seed ){
    this._hash = Cell.genarateHash( seed.x, seed.y )
    this._polygon= new ConvexPolygon( contourPoints );
    this._seed = seed;
    this._isValid = (this._polygon.isValid() && !!this._hash);
  }


  /**
  * Return if this cell is valid
  * @return {Boolean} true if valid, false if not
  */
  isValid(){
    return this._isValid;
  }


  /**
  * Get the hash of this cell
  * @return {String}
  */
  getHash(){
    return this._hash;
  }


  /**
  * Get the seed of the cell.
  * @return {Object} should be of form {x: Number, y: Number, seedIndex: Number}
  */
  getSeed(){
    return this._seed;
  }


  /**
  * Get the polygon of this cell
  * @return {ConvexPolygon}
  */
  getPolygon(){
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
  static genarateHash( n1, n2 ){
    //return md5( new Float32Array([f1, f2]) );
    var base = 16;
    //var hash = (i1 < base ? "0" : '') + i1.toString(16) +
    //           (i2 < base ? "0" : '') + i2.toString(16)

    //
    var hash = n1.toString() + "_" + n2.toString();
    return hash;
  }


  /**
  * Get the area of the cell (calls the polygon method)
  * @return {Number} the area
  */
  getArea(){
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
  intersectWithCell( anotherCell ){
    return this._polygon.getCellReplacementIntersection( anotherCell.getPolygon() );
  }


  /**
  * Compare if this cell has the same polygon as another cell.
  * Read the doc of ConvexPolygon.getPolygon() for more info.
  * @param {Cell} anotherCell - another cell to compare the polygon with.
  * @return {Boolean} true if the same polygon, false if not
  */
  hasSamePolygon( otherCell ){
    return this._polygon.isSame( otherCell.getPolygon() );
  }

}

export { Cell };
