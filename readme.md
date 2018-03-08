# 2D natural neighbor interpolation (nni)
For Node and browsers.  

Some info [about nni](http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-natural-neighbor-works.htm).  
The *nii* method is used to interpolate sparse data point along a regular 2D grid surface. This methods work with Voronoi cells and area stilling and is said to give results that are close to the "real solution". The only issue with this method is that it's very computational, and thus, quite slow.

## Concept
Sparse point interpolations need:  
- a set of sparse points. Here, we will call them *seeds*
- an output size (width and height)

This implementation performs the interpolation in two steps:
- build the interpolation map for each output pixel, aka the weight of each seeds involved in each output pixel
- use the interpolation map to create the output image

The purpose of creating the map at an intermediary step is to be able to reuse it for the same configuration but with different seed values, because what takes a while to generate is this interpolation map.

## Basic usage
```Javascript
// a list of seeds
var seeds = [
  {x: 100, y: 100,  value: Math.random()*255 },
  {x: 300, y: 100,  value: Math.random()*255  },
  {x: 300, y: 300, value: Math.random()*255 },
  {x: 100, y: 300, value: Math.random()*255  },
  {x: 200 , y: 200 , value: Math.random()*255 }
]

// creating the nni interpolator instance
var nnInter = new natninter.Interpolator();

// setting the output size
nnInter.setOutputSize(256, 256);

// add a list of seeds
nnInter.addSeeds( seeds );

// Create the interpolation map
// (this may take some seconds, start with a small image to benchmark it)
nnInter.generateMap();

// generate the output image witha  nice interpolation
var output = nnInter.generateImage();
```

The `output` is actually an object of the form:
```Javascript
{
  _data: Float32Array,
  _metadata: {
    width: Number,
    heigth: Number
  }
}
```
Where `_data` is a 1D array of length `height`x`width` arrange in a row major way.

## Examples
See the folder `examples` for a Node and a browser example.

## Generate a map with and save it locally
`natninter` comes with a *node* executable to generate a samplig-map file in JSON format. In a terminal, find the executable `generatemap.js`, located in the `bin` folder. Then:

```bash
./generatemap.js --seeds=your/own/seedsfile.json -width=400 --height=400 --output=somehwere/map.json
```
You can use the seed file in `bin/samples/seeds1.json` and use it as a a model to write your own seed file.


# Special thanks
Since I didn't want to reinvent the wheel and that I trust some good developers out there, here are the dependencies *natninter* uses:
- [Compute Voronoi diagram](https://github.com/gorhill/Javascript-Voronoi) - damn fast!
- [Polygon area computation](https://github.com/math-utils/area-polygon) - super simple a pretty fast too!
