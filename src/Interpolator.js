import Voronoi from 'voronoi';
import { CellCollection } from './CellCollection.js';



class Interpolator {
  constructor(){
    this._seeds = [];
    this._output = { width: 256, height: 256 };
    this._samplingMap = [];
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


  hasAllSeedsInside(){
    var size = this._output;
    return this._seeds.every( function(s){
      return ( s.x>=0 && s.y>=0 && s.x<size.width && s.y<size.width );
    })
  }


  /**
  * Compute the sampling map. Automatically called by the update method when the
  * map was never computed or when a seed have been added since.
  * Though this method is not private and can be called to force recomputing the map
  */
  computeMap(){
    if( !this.hasAllSeedsInside() ){
      console.log( 'ERR: some seeds are outside of the image. Use .setOutputSize() to change the output size or modify the seed.' );
      return;
    }
    this._generateSeedCells();
    console.log( this._seedCellCollection );

    var w = this._output.width;
    var h = this._output.height;
    this._samplingMap = new Array( w*h );

    // for each pixel of the output image
    for(var i=0; i<w; i++){
      for(var j=0; j<h; j++){

        var seedIndex = this.isAtSeedPosition(i, j);
        var index1D = i + j*w;

        if( seedIndex === -1 ){
          //this._samplingMap[ index1D ] = this._generatePixelCells(i, j);
          var pixelCellCollection = this._generatePixelCells(i, j);
          var stolenAreaData = this._seedCellCollection.getStolenAreaInfo( pixelCellCollection )
          this._samplingMap[ index1D ] = stolenAreaData;
        }else{
          this._samplingMap[ index1D ] = [{seedIndex: seedIndex, stolenRatio: 1}];
        }
      }
    }


    console.log( this._samplingMap );
    // TODO: finish map for seed positions
  }


  getMap( m ){
    return this._samplingMap;
  }


  setMap( map ){
    if( map.length === (this._output.width * this._output.height)){
      this._samplingMap = map;
    }else{
      console.log("The sampling map must be an 1D array of size output.x*output.y");
    }
  }

  /**
  * is the given position at the position of a seed?
  * @param {Number} i - position along x axis
  * @param {Number} j - position along y axis
  * @return {Boolean} -1 if not, of the index of the seed if yes
  */
  isAtSeedPosition(i, j){
    return this._seeds.findIndex(function(elem){ return (elem.x==i && elem.y==j)})
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
    var seedVoronoiDiagram = voronoi.compute(sites, bbox);

    this._seedCellCollection = new CellCollection();
    this._seedCellCollection.buildFromVoronoiDiagram( seedVoronoiDiagram );
  }


  _generatePixelCells(i, j){
    var that = this;
    var voronoi = new Voronoi();
    var bbox = {xl: 0, xr: this._output.width, yt: 0, yb: this._output.height};
    var sites = this._seeds.map( function( s, i ){ return {x: s.x, y: s.y, seedIndex: i} });
    sites.push( {x: i, y: j, seedIndex: -1} )
    var pixelVoronoiDiagram = voronoi.compute(sites, bbox);

    var pixelCellCollection = new CellCollection();
    pixelCellCollection.buildFromVoronoiDiagram( pixelVoronoiDiagram );

    if( pixelCellCollection._cellArray.length !== (this._seeds.length + 1))
      console.log("ERR");

    return pixelCellCollection;
  }


  /**
  * Generate the output image as a floating points 1D array representing a 2D (1band)
  * image.
  *
  */
  generate(){
    var l = this._output.width * this._output.height;
    var outImg = new Float32Array( l );
    var seeds = this._seeds;
    var map = this._samplingMap;

    for(var i=0; i<l; i++){
      var pixelMap = map[i];
      var sum = 0;

      for(var m=0; m<pixelMap.length; m++){
        sum += (pixelMap[m].stolenRatio *  seeds[ pixelMap[m].seedIndex ].value  )
      }

      outImg[i] = sum;
    }

    return {
      _metadata: {
        width: this._output.width,
        height: this._output.height
      },
      _data: outImg
    }
  }


}

export { Interpolator };
