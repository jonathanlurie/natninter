import { Cell } from './Cell.js';

class CellCollection {

  constructor(){
    this._cells = {}
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
      }
    }
  }


}


export { CellCollection };
