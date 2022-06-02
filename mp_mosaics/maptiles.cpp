/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    vector<Point<3>> tile_avg_color_points;

    // a constant time lookup dictionary for retrieving the correct tile from its average color
    map<Point<3>, TileImage*> avg_color_to_tile;

    tile_avg_color_points.reserve(theTiles.size());
    for (size_t i = 0; i < theTiles.size(); ++i) {
        Point<3> curr_tile_avg_color_point = convertToXYZ(theTiles[i].getAverageColor());
        tile_avg_color_points.push_back(curr_tile_avg_color_point);
        avg_color_to_tile.insert({curr_tile_avg_color_point, &theTiles[i]});
    }

    // construct a kd_tree (k = 3) which is used for determining the nearest neighbor 
    KDTree<3> tree(tile_avg_color_points);

    int canvas_rows = theSource.getRows();
    int canvas_cols = theSource.getColumns();
    MosaicCanvas* canvas = new MosaicCanvas(canvas_rows, canvas_cols);
     
    for (int j = 0; j < canvas_rows; ++j) {
        for (int i = 0; i < canvas_cols; ++i) {
            Point<3> source_region_avg_color_point = convertToXYZ(theSource.getRegionColor(i, j));
            Point<3> nearest_neighbor_avg_color = tree.findNearestNeighbor(source_region_avg_color_point);
            // map the nearest neighbor average color to its tile
            canvas->setTile(i, j, avg_color_to_tile[nearest_neighbor_avg_color]);
        }
    }

    return canvas;
}

