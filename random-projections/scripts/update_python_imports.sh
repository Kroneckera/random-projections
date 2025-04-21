#!/bin/bash

# Script to update Python imports from polygon_projection to projection

# Process all Python files
for file in $(find /Users/azimin/Programming/random_proj/random-projections/python -name "*.py"); do
    echo "Processing $file..."
    # Update "from polygon_projection" to "from projection"
    sed -i '' 's/from polygon_projection/from projection/g' "$file"
    # Update "import polygon_projection" to "import projection"
    sed -i '' 's/import polygon_projection/import projection/g' "$file"
    # Update references to polygon_projection.XXX to projection.XXX
    sed -i '' 's/polygon_projection\./projection\./g' "$file"
done

# Process all test Python files
for file in $(find /Users/azimin/Programming/random_proj/random-projections/python/tests -name "*.py"); do
    echo "Processing test file $file..."
    # Update "from polygon_projection" to "from projection"
    sed -i '' 's/from polygon_projection/from projection/g' "$file"
    # Update "import polygon_projection" to "import projection"
    sed -i '' 's/import polygon_projection/import projection/g' "$file"
    # Update references to polygon_projection.XXX to projection.XXX
    sed -i '' 's/polygon_projection\./projection\./g' "$file"
done

echo "Python imports updated successfully!"