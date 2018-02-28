import Voronoi from 'voronoi';

class Interpolator {
  constructor(){
    this._seeds = [];
    this._output = { width: 256, height: 256 };
    this._weightMap = [];
    this._recomputeMap = true;

    // voronoi diagram of seeds only
    this._seedVoronoiDiagram = null;
  }


  /**
  * Removes all the seeds
  */
  cleanSeeds(){
    this._seeds = [];
    this._recomputeMap = true;
  }

  /**
  * Add a seed to the interpolator. Note that the seed is not copied, only its reference is.
  * This is convenient for updating the value of the seed from the outside and ecompute
  * the interpolation without having to recompute the whole weight map.
  * @param {Object} seed - of form {x: Number, y: Number, value: Number}
  */
  addSeed( seed ){
    if( "x" in seed && "y" in seed && "value" in seed ){
      this._seeds.push( seed );
      this._recomputeMap = true;
    }
  }


  /**
  * Add an array of seed
  * @param {Array} seedArr - array of seeds, where each seed is of form {x: Number, y: Number, value: Number}
  */
  addSeeds( seedArr ){
    for(var i=0; i<seedArr.length; i++){
      this.addSeed( seedArr[i] );
    }
  }


  /**
  * Define the size of the output image
  */
  setOutputSize( width, height ){
    if( width > 0 && height > 0 ){
      this._output.width = width;
      this._output.height = height;
    }
  }


  /**
  * Compute the sampling map. Automatically called by the update method when the
  * map was never computed or when a seed have been added since.
  * Though this method is not private and can be called to force recomputing the map
  */
  computeMap(){
    this._generateSeedVoronoiDiagram();
  }


  /**
  * [PRIVATE]
  * Generate the voronoi diagram where sites are only seeds
  */
  _generateSeedVoronoiDiagram(){
    var voronoi = new Voronoi();
    var bbox = {xl: 0, xr: this._output.width, yt: 0, yb: this._output.height};

    var sites = new Array( this._seeds.length );
    for(var i=0; i<sites.length; i++){
      sites[i] = {x: this._seeds[i].x, y: this._seeds[i].y, seedIndex: i}
    }

    // we could copy the ref but i am not sure what Voronoi does with it.
    // I must check before messing with ref.
    // var sites = this._seeds

    // a 'vertex' is an object exhibiting 'x' and 'y' properties. The
    // Voronoi object will add a unique 'voronoiId' property to all
    // sites. The 'voronoiId' can be used as a key to lookup the associated cell
    // in diagram.cells.

    console.log( sites );
    this._seedVoronoiDiagram = voronoi.compute(sites, bbox);
    console.log( this._seedVoronoiDiagram );
  }



}

export { Interpolator };
