import Voronoi from 'voronoi';
import { CellCollection } from './CellCollection.js';



class Interpolator {
  constructor(){
    this._seeds = [];
    this._output = { width: 256, height: 256 };
    this._weightMap = [];
    this._recomputeMap = true;

    // voronoi diagram of seeds only
    this._seedCellCollection = null;
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
    //console.time("0")
    //for(var i=0; i<512*512; i++)
    this._generateSeedCells();
    //console.timeEnd("0")
  }


  /**
  * [PRIVATE]
  * Generate the voronoi diagram where sites are only seeds
  */
  _generateSeedCells(){
    var that = this;
    var voronoi = new Voronoi();
    var bbox = {xl: 0, xr: this._output.width, yt: 0, yb: this._output.height};
    var sites = this._seeds.map( function( s, i ){ return {x: s.x, y: s.y, seedIndex: i} });

    // we could copy the ref but i am not sure what Voronoi does with it.
    // I must check before messing with ref.
    // var sites = this._seeds

    // a 'vertex' is an object exhibiting 'x' and 'y' properties. The
    // Voronoi object will add a unique 'voronoiId' property to all
    // sites. The 'voronoiId' can be used as a key to lookup the associated cell
    // in diagram.cells.

    var seedVoronoiDiagram = voronoi.compute(sites, bbox);

    this._seedCellCollection = new CellCollection();
    this._seedCellCollection.buildFromVoronoiDiagram( seedVoronoiDiagram );
  }



}

export { Interpolator };
