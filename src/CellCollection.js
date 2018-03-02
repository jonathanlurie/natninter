import { Cell } from './Cell.js';

class CellCollection {

  constructor(){
    this._cells = {}
    this._cellArray = [];
  }

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
  */
  getCell( hash ){
    if(hash in this._cells)
      return this._cells[hash];
  }


  getUnindexedCells(){
    var  unindexedCells = [];

    for(var i=0; i<this._cellArray.length; i++){
      var cell = this._cellArray[ i ]
      if( cell.getSeed().seedIndex === -1 )
        unindexedCells.push( cell );
    }
    /*
    for( var hash in this._cells ){
      var cell = this._cells[ hash ]
      if( cell.getSeed().seedIndex === -1 )
        unindexedCells.push( cell );
    }
    */
    return unindexedCells;
  }

  /**
  * Get the list of cells that are in the cell given in arg and that are NOT in
  * _this_ collection.
  * @param {CellCollection} anotherCellCollection - another cell collection, that presumably has cells that _this_ collection does not.
  */
  getAdditionalCells( anotherCellCollection ){
    var anotherCellHashes = anotherCellCollection.getCellHashes()
    var extraCells = [];

    for(var i=0; i<anotherCellHashes.length; i++){
      if( !(anotherCellHashes[i] in this._cells) ){
        extraCells.push( anotherCellCollection.getCell(anotherCellHashes[i]) );
      }
    }
    return extraCells;
  }


  getStolenAreaInfo( anotherCellCollection ){
    // this is the added cell from the original
    var addedCells = anotherCellCollection.getUnindexedCells();

    if(addedCells.length === 0)
      return null;

    var addedCell = addedCells[0];
    var addedCellArea = addedCell.getArea();

    var modifiedCells = [];

    for( var hash in this._cells ){
      if( hash === addedCell.getHash() )
        continue;

      var originalCell = this._cells[ hash ];
      var newCell = anotherCellCollection.getCell( hash );
      var isSame = originalCell.hasSamePolygon( newCell );

      if( !isSame ){
        var stolenPoly = addedCell.intersectWithCell( originalCell );
        var stolenArea = stolenPoly.getArea();
        var stolenRatio = stolenArea / addedCellArea;
        modifiedCells.push( {cell: originalCell, stolenRatio: stolenRatio} );
      }
    }

    console.log( modifiedCells );
    // TODO: test that!!!
    return modifiedCells;


  }


}


export { CellCollection };
