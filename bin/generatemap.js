#!/usr/bin/env node

var fs = require("fs");
var ArgParser = require("./tools/ArgParser");
var argParser = new ArgParser();

const natninter = require("../dist/natninter.js");

function help(){
  var helpStr = `
  ______________________________________________________________________________

  The executable ${argParser.getAppName()} creates a reasampling map file using
  the algorithm of the Natural Neighbor Interpolation. The ressampling map is a
  JSON file that can then be used by other appilication usinf 'natninter'.
  The following arguments are required:

    --seeds=/path/to/seedfile/json   - The path to a seed description file
    --width=256                      - The width of the resampling map
    --heigth=256                     - The width of the resampling map
    --output=/output/map.json        - The path to the output resampling map

  Alternatively:
    --help                           - Print this help menu

  ______________________________________________________________________________

  `;
  console.log( helpStr );
  process.exit();
}


function printPercentBar( percent ){
  process.stdout.clearLine();  // clear current text
  process.stdout.cursorTo(0);

  var status = new Array(Math.round(percent/5)+1).join("▓") +
               new Array(20 - Math.round(percent/5) + 1).join("░") + ' ' + percent+"%";
  process.stdout.write(status);
}


function printStatusUpdate( message ){
  process.stdout.clearLine();  // clear current text
  process.stdout.cursorTo(0);
  process.stdout.write("📣 " + message + "\n");
}




// print help menu if no options
if( argParser.getNumberOfArgs() === 0 ){
  help();
}


// do we have to print the help?
try{
  var helpVal = argParser.getArgValue("help");
  if( helpVal )
    help();
}catch(e){}




var seedFile = null;
var width = null;
var height = null;
var output = null;


try{
  seedFile = argParser.getArgValue("seeds");
  width = argParser.getArgValue("width");
  height = argParser.getArgValue("height");
  output = argParser.getArgValue("output");
}catch( e ){
  printStatusUpdate( '\n\tERROR: '+  e  )
  help();
}

process.stdout.write("\n");

// getting the seed from the seed file
var seeds = null;
try{
  seeds = JSON.parse( fs.readFileSync( seedFile ) );
}catch(e){
  printStatusUpdate( e.message );
  process.exit();
}


var nnInter = new natninter.Interpolator();

nnInter.onProgress( function( percent, message){
  printPercentBar( percent );
})

nnInter.setOutputSize(width, height);
nnInter.addSeeds( seeds );

nnInter.generateMap();
printStatusUpdate("Map generated.\n")


var map = null;

try{
  map = JSON.stringify( nnInter.getMap() );
}catch(e){
  printStatusUpdate( e.message );
  process.exit();
}


try{
  fs.writeFileSync(output, map);
  printStatusUpdate( "Map file exported in " + output );
}catch(e){
  printStatusUpdate( e.message );
  process.exit();
}

process.exit();
