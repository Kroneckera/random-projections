#!/bin/bash

# Script to update C++ include directives from polygon_projection/file.h to file.h

# Process all C++ source files
for file in $(find /Users/azimin/Programming/random_proj/random-projections/cpp -name "*.cpp"); do
    echo "Processing $file..."
    # Update #include "polygon_projection/file.h" to #include "file.h"
    sed -i '' 's/#include "polygon_projection\//#include "/g' "$file"
    # Update include guardss from POLYGON_PROJECTION_FILE_H to PROJECTION_FILE_H
    sed -i '' 's/POLYGON_PROJECTION_/PROJECTION_/g' "$file"
done

# Process all C++ header files
for file in $(find /Users/azimin/Programming/random_proj/random-projections/cpp -name "*.h"); do
    echo "Processing $file..."
    # Update #include "polygon_projection/file.h" to #include "file.h"
    sed -i '' 's/#include "polygon_projection\//#include "/g' "$file"
    # Update include guards from POLYGON_PROJECTION_FILE_H to PROJECTION_FILE_H
    sed -i '' 's/POLYGON_PROJECTION_/PROJECTION_/g' "$file"
done

# Update namespace from polygon_projection to projection
for file in $(find /Users/azimin/Programming/random_proj/random-projections/cpp -name "*.cpp" -o -name "*.h"); do
    echo "Updating namespace in $file..."
    sed -i '' 's/namespace polygon_projection/namespace projection/g' "$file"
    sed -i '' 's/polygon_projection::/projection::/g' "$file"
done

echo "C++ includes updated successfully!"