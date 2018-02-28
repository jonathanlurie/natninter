**WIP**  

# 2D natural neighbor interpolation (nni)
Some info [here](http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-natural-neighbor-works.htm).  
The *nii* method is used to interpolate sparse data point along a regular 2D grid surface. This methods work with Voronoi cells and area stilling.

Since I didn't want to reinvent the wheel and that I trust some good developers out there, here are the dependencies *natninter* uses:
- [Compute Voronoi diagram](https://github.com/gorhill/Javascript-Voronoi)
- [Polygon clipping](https://github.com/w8r/GreinerHormann)
- [Polygon area computation](https://github.com/math-utils/area-polygon)
- [Compute convex hull](https://github.com/mikolalysenko/convex-hull)
