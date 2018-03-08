
/**
 * VectorTools is not instanciable and provide some static functions for
 * computing things about vectors.
 */
class VectorTools {

  /**
   * performs a cross product between v1 and v2.
   * @param  {Array} v1 - a vector [x, y, z]. To use with 2D vectors, just use [x, y, 0]
   * @param  {Array} v2 - a vector [x, y, z]. To use with 2D vectors, just use [x, y, 0]
   * @param  {Boolean} normalize - will force normalization of the output vector if true (default: false)
   * @return {Array} a vector [x, y, z], result of the cross product
   */
  static crossProduct(v1, v2, normalize=false){
    var normalVector = [
        v1[1] * v2[2] - v1[2] * v2[1],
        (v1[0] * v2[2] - v1[2] * v2[0] ) * (-1),
        v1[0] * v2[1] - v1[1] * v2[0]
    ];

    if(normalize)
        normalVector = VectorTools.normalize(normalVector);

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
  static dotProduct(v1, v2){
    if(v1.length != v2.length){
      console.log("ERROR: v1 and v2 must be the same size to compute a dot product");
      return null;
    }

    var sum = 0;

    for(var i=0; i<v1.length; i++){
      sum += (v1[i]*v2[i]);
    }

    return sum;
  }


 /**
  * Normalizes a 3D vector
  * @param  {Array} v - 3D vector to normalize
  * @return {Array} the normalized 3D vector
  */
  static normalize(v){
    var n = VectorTools.getNorm(v);
    var normalizedV = [v[0]/n, v[1]/n, v[2]/n];
    return normalizedV;
  }


 /**
  * return the norm (length) of a vector [x, y, z]
  * @param  {Array} v - 3D vector to get the norm of
  * @return {Number} the norm
  */
  static getNorm(v){
    return Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
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
  static rotate(v, m){
    var vRot = [
      v[0]*m[0][0] + v[1]*m[0][1] + v[2]*m[0][2],
      v[0]*m[1][0] + v[1]*m[1][1] + v[2]*m[1][2],
      v[0]*m[2][0] + v[1]*m[2][1] + v[2]*m[2][2],
    ];

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
  static getAnglePoints(p1, p2, p3){

    //  the configuration is like that:
    //
    //  p1-----p2
    //        /
    //       /
    //      p3

    var v_p2p1 = [
      p1[0] - p2[0],
      p1[1] - p2[1],
      p1[2] - p2[2],
    ];

    var v_p2p3 = [
      p3[0] - p2[0],
      p3[1] - p2[1],
      p3[2] - p2[2],
    ];

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
  static vector2DCrossing(u1, u2, v1, v2){
    var meet = null;
    var s1x = u2[0] - u1[0];
    var s1y = u2[1] - u1[1];
    var s2x = v2[0] - v1[0];
    var s2y = v2[1] - v1[1];

    var s = (-s1y * (u1[0] - v1[0]) + s1x * (u1[1] - v1[1])) / (-s2x * s1y + s1x * s2y);
    var t = ( s2x * (u1[1] - v1[1]) - s2y * (u1[0] - v1[0])) / (-s2x * s1y + s1x * s2y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
      meet = [
        u1[0] + (t * s1x),
        u1[1] + (t * s1y)
      ]
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
  static pointToSegmentDistance(p, u1, u2) {
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
    }
    else if (param > 1) {
      xx = u2[0];
      yy = u2[1];
    }
    else {
      xx = u1[0] + param * C;
      yy = u1[1] + param * D;
    }

    var dx = p[0] - xx;
    var dy = p[1] - yy;
    return Math.sqrt(dx * dx + dy * dy);
  }


}

export { VectorTools };
