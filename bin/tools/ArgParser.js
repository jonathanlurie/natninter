const path = require('path');

class ArgParser {
  constructor(){
    //console.log( process.argv );

    this._parsedArgs = {}
    this._appName = path.basename( process.argv[1] );
    this._parseArgs();
  }


  _parseArgs(){
    // reinit the thing
    this._parsedArgs = {};
    var rgx = /([-]+)([a-zA-Z]+)(=(\S+))?/;

    var argz = process.argv

    for(var i=0; i<argz.length; i++){
        var arg = argz[i];
        var match = rgx.exec(arg);

        if( match ){
          var argName = match[2]; // group 2
          var argValue = ArgParser._castValue( match[4] )  // group 4
          this._parsedArgs[ argName ] = argValue;
          //console.log( argName);
        }

    }
  }


  static _castValue( value ){
    if( value === undefined ){
      return true;
    }

    if( value === "false" || value === "False" || value === "FALSE"){
      return false;
    }

    if( value === "true" || value === "True" || value === "TRUE"){
      return true;
    }

    var numberVal = parseFloat( value, 10 );

    if( isNaN(numberVal) ){
      return value;
    }else{
      return numberVal;
    }
  }


  getArgValue( argName ){
    if( argName in this._parsedArgs){
      return this._parsedArgs[argName];
    }else{
      throw "The argument --" + argName + " is missing."
    }
  }


  getAppName(){
    return this._appName;
  }


  getNumberOfArgs(){
    return Object.keys( this._parsedArgs ).length;
  }
}

module.exports = ArgParser;
