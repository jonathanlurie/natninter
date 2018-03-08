import { Cell } from './Cell.js';


/**
* A CellCollection stores Cell objects. It is the interface between the Voronoi diagram
* and the Cells/Polygons.
*/
class CellCollection {

  /**
  * The constructor does not take any param.
  */
  constructor(){
    this._cells = {}
    this._cellArray = [];
    this._addedSeedHash = null;
  }


  /**
  * Builds the collection of cell using a voronoi diagram object from the Javascript-Voronoi
  * library by Gorhill ( http://bit.ly/2FoapPI )
  */
  buildFromVoronoiDiagram( vd ){
    var vCells = vd.cells;

    // for each cell
    for(var i=0; i<vCells.length; i++){
      var cellhes = vCells[i].halfedges;
      var points = [];

      for(var j=0; j<cellhes.length; j++){
        points.push( cellhes[j].edge.va );
        points.push( cellhes[j].edge.vb );
      }

      // build a cell
      var cell = new Cell( points, vCells[i].site );
      if( cell.isValid() ){
        this._cells[ cell.getHash() ] = cell;
        this._cellArray.push( cell );
      }
    }
  }


  /**
  * Get the list of hashes (ID) of cells in the collection
  * @return {Array} all the hashes
  */
  getCellHashes(){
    return Object.keys( this._cells );
  }


  /**
  * Get the cell of _this_ collection that has this hash
  * @param {String} hash - the unique hash of the cell
  * @return {Cell} a Cell instance of null if hash not found
  */
  getCell( hash ){
    if(hash in this._cells)
      return this._cells[hash];
    else
      return null;
  }


  /**
  * Store the hash of a given position of a seed, meaning, store a reference to a cell.
  * This seed is the one added to the 'another' diagram, when we introduce a pixel as a seed.
  * For instance, the seed-only cell collection does not have a reference to a seed.
  * @param {Number} x - the x position of the seed to reference
  * @param {Number} y - the y position of the seed to reference
  */
  referenceSeed(x, y){
    this._addedSeedHash = Cell.genarateHash(x, y);
  }


  /**
  * Get the cell that was referenced by referenceSeed()
  * @return {Cell}
  */
  getReferencedSeedCell(){
    return this._cells[ this._addedSeedHash ];
  }


  /**
  * When comparing a seed-only CellCollection with a seed-and-one-pixel CellCollection
  * with getStolenAreaInfo(), we get a list of original cell index (from the seed-only collection)
  * as well as the area ratio that comes from them to build the pixel-based cell
  * @param {CellCollection} anotherCellCollection - a cell collection with the original seeds + 1 pixel as a seed
  * @return {Array} list of {seedIndex: Number, weight: Number}
  */
  getStolenAreaInfo( anotherCellCollection ){
    var addedCell = anotherCellCollection.getReferencedSeedCell();
    var addedCellArea = addedCell.getArea();
    var modifiedCells = [];
    var stolenTotal = 0;

    for( var hash in this._cells ){
      if( hash === addedCell.getHash() )
        continue;

      var originalCell = this._cells[ hash ];
      var newCell = anotherCellCollection.getCell( hash );
      var isSame = originalCell.hasSamePolygon( newCell );

      if( !isSame ){
        var stolenRatio = 0;
        var stolenPoly = addedCell.intersectWithCell( originalCell );

        if( stolenPoly.isValid() ){
          var stolenArea = stolenPoly.getArea();
          stolenRatio = stolenArea / addedCell.getArea();
          stolenTotal += stolenRatio;
        }

        stolenRatio = Math.round( stolenRatio * 10000) / 10000;

        //modifiedCells.push( {cell: originalCell, stolenRatio: stolenRatio} );
        modifiedCells.push( {seedIndex: originalCell.getSeed().seedIndex, weight: stolenRatio} );
      }
    }

    return modifiedCells;
  }

}


export { CellCollection };
