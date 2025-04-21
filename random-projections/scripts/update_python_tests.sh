#!/bin/bash

# Script to update Python test imports from polygon_projection to projection

# Process all Python test files
for file in $(find /Users/azimin/Programming/random_proj/random-projections/python/tests -name "*.py"); do
    echo "Processing test file $file..."
    
    # Update test script header comment
    sed -i '' 's/the polygon_projection module/the projection module/g' "$file"
    sed -i '' 's/for polygon-projection/for projection/g' "$file"
    
    # Update environment variable name for forcing fallback
    sed -i '' 's/POLYGON_PROJECTION_FORCE_FALLBACK/PROJECTION_FORCE_FALLBACK/g' "$file"
    
    # Update import statements
    sed -i '' 's/polygon_projection as pp/projection as pp/g' "$file"
    sed -i '' 's/from polygon_projection/from projection/g' "$file"
    sed -i '' 's/import polygon_projection/import projection/g' "$file"
done

echo "Python test imports updated successfully!"