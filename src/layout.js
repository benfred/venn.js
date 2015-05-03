(function(venn) {
    /** given a list of set objects, and their corresponding overlaps.
    updates the (x, y, radius) attribute on each set such that their positions
    roughly correspond to the desired overlaps */
    venn.venn = function(areas, parameters) {
        parameters = parameters || {};
        parameters.maxIterations = parameters.maxIterations || 500;
        var lossFunction = parameters.lossFunction || venn.lossFunction;
        var initialLayout = parameters.initialLayout || venn.greedyLayout;
        var fmin = parameters.fmin || venn.fmin;

        // initial layout is done greedily
        var circles = initialLayout(areas);

        // transform x/y coordinates to a vector to optimize
        var initial = [], setids = [], setid;
        for (setid in circles) {
            if (circles.hasOwnProperty(setid)) {
                initial.push(circles[setid].x);
                initial.push(circles[setid].y);
                setids.push(setid);
            }
        }

        // optimize initial layout from our loss function
        var totalFunctionCalls = 0;
        var solution = fmin(
            function(values) {
                totalFunctionCalls += 1;
                var current = {};
                for (var i = 0; i < setids.length; ++i) {
                    var setid = setids[i];
                    current[setid] = {x: values[2 * i],
                                      y: values[2 * i + 1],
                                      radius : circles[setid].radius,
                                     // size : circles[setid].size
                                     };
                }
                return lossFunction(current, areas);
            },
            initial,
            parameters);

        // transform solution vector back to x/y points
        var positions = solution.solution;
        for (var i = 0; i < setids.length; ++i) {
            setid = setids[i];
            circles[setid].x = positions[2 * i];
            circles[setid].y = positions[2 * i + 1];
        }

        return circles;
    };
    
    var SMALL = 1e-10;

    /** Returns the distance necessary for two circles of radius r1 + r2 to
    have the overlap area 'overlap' */
    venn.distanceFromIntersectArea = function(r1, r2, overlap) {
        // handle complete overlapped circles
        if (Math.min(r1, r2) * Math.min(r1,r2) * Math.PI <= overlap + SMALL) {
            return Math.abs(r1 - r2);
        }

        return venn.bisect(function(distance) {
            return venn.circleOverlap(r1, r2, distance) - overlap;
        }, 0, r1 + r2);
    };

    /// gets a matrix of euclidean distances between all sets in venn diagram
    venn.getDistanceMatrix = function(areas, sets, setids) {
        // initialize an empty distance matrix between all the points
        var distances = [];
        for (var i = 0; i < sets.length; ++i) {
            distances.push([]);
            for (var j = 0; j < sets.length; ++j) {
                distances[i].push(0);
            }
        }

        // compute distances between all the points
        for (i = 0; i < areas.length; ++i) {
            var current = areas[i];
            if (current.sets.length !== 2) {
                continue;
            }

            var left = setids[current.sets[0]],
                right = setids[current.sets[1]],
                r1 = Math.sqrt(sets[left].size / Math.PI),
                r2 = Math.sqrt(sets[right].size / Math.PI),
                distance = venn.distanceFromIntersectArea(r1, r2, current.size);

            distances[left][right] = distances[right][left] = distance;
        }
        return distances;
    };

    /** Lays out a Venn diagram greedily, going from most overlapped sets to
    least overlapped, attempting to position each new set such that the
    overlapping areas to already positioned sets are basically right */
    venn.greedyLayout = function(areas) {
        // define a circle for each set
        var circles = {}, setOverlaps = {}, set;
        for (var i = 0; i < areas.length; ++i) {
            var area = areas[i];
            if (area.sets.length == 1) {
                set = area.sets[0];
                circles[set] = {x: 1e10, y: 1e10,
                                rowid: circles.length,
                                size: area.size,
                                radius: Math.sqrt(area.size / Math.PI)};
                setOverlaps[set] = [];
            }
        }
        areas = areas.filter(function(a) { return a.sets.length == 2; });

        // map each set to a list of all the other sets that overlap it
        for (i = 0; i < areas.length; ++i) {
            var current = areas[i];
            var weight = current.hasOwnProperty('weight') ? current.weight : 1.0;
            var left = current.sets[0], right = current.sets[1];

            // completely overlapped circles shouldn't be positioned early here
            if (current.size + SMALL >= Math.min(circles[left].size,
                                                 circles[right].size)) {
                weight = 0;
            }

            setOverlaps[left].push ({set:right, size:current.size, weight:weight});
            setOverlaps[right].push({set:left,  size:current.size, weight:weight});
        }

        // get list of most overlapped sets
        var mostOverlapped = [];
        for (set in setOverlaps) {
            if (setOverlaps.hasOwnProperty(set)) {
                var size = 0;
                for (i = 0; i < setOverlaps[set].length; ++i) {
                    size += setOverlaps[set][i].size * setOverlaps[set][i].weight;
                }

                mostOverlapped.push({set: set, size:size});
            }
        }

        // sort by size desc
        function sortOrder(a,b) {
            return b.size - a.size;
        }
        mostOverlapped.sort(sortOrder);

        // keep track of what sets have been laid out
        var positioned = {};
        function isPositioned(element) {
            return element.set in positioned;
        }

        // adds a point to the output
        function positionSet(point, index) {
            circles[index].x = point.x;
            circles[index].y = point.y;
            positioned[index] = true;
        }

        // add most overlapped set at (0,0)
        positionSet({x: 0, y: 0}, mostOverlapped[0].set);

        // get distances between all points. TODO, necessary?
        // answer: probably not
        // var distances = venn.getDistanceMatrix(circles, areas);
        for (i = 1; i < mostOverlapped.length; ++i) {
            var setIndex = mostOverlapped[i].set,
                overlap = setOverlaps[setIndex].filter(isPositioned);
            set = circles[setIndex];
            overlap.sort(sortOrder);

            if (overlap.length === 0) {
                throw "Need overlap information for set " + JSON.stringify( set );
            }

            var points = [];
            for (var j = 0; j < overlap.length; ++j) {
                // get appropriate distance from most overlapped already added set
                var p1 = circles[overlap[j].set],
                    d1 = venn.distanceFromIntersectArea(set.radius, p1.radius,
                                                        overlap[j].size);

                // sample positions at 90 degrees for maximum aesthetics
                points.push({x : p1.x + d1, y : p1.y});
                points.push({x : p1.x - d1, y : p1.y});
                points.push({y : p1.y + d1, x : p1.x});
                points.push({y : p1.y - d1, x : p1.x});

                // if we have at least 2 overlaps, then figure out where the
                // set should be positioned analytically and try those too
                for (var k = j + 1; k < overlap.length; ++k) {
                    var p2 = circles[overlap[k].set],
                        d2 = venn.distanceFromIntersectArea(set.radius, p2.radius,
                                                            overlap[k].size);

                    var extraPoints = venn.circleCircleIntersection(
                        { x: p1.x, y: p1.y, radius: d1},
                        { x: p2.x, y: p2.y, radius: d2});

                    for (var l = 0; l < extraPoints.length; ++l) {
                        points.push(extraPoints[l]);
                    }
                }
            }

            // we have some candidate positions for the set, examine loss
            // at each position to figure out where to put it at
            var bestLoss = 1e50, bestPoint = points[0];
            for (j = 0; j < points.length; ++j) {
                circles[setIndex].x = points[j].x;
                circles[setIndex].y = points[j].y;
                var loss = venn.lossFunction(circles, areas);
                if (loss < bestLoss) {
                    bestLoss = loss;
                    bestPoint = points[j];
                }
            }

            positionSet(bestPoint, setIndex);
        }

        return circles;
    };

    /// Uses multidimensional scaling to approximate a first layout here
    venn.classicMDSLayout = function(areas) {
        // bidirectionally map sets to a rowid  (so we can create a matrix)
        var sets = [], setids = {};
        for (var i = 0; i < areas.length; ++i ) {
            var area = areas[i];
            if (area.sets.length == 1) {
                setids[area.sets[0]] = sets.length;
                sets.push(area);
            }
        }

        // get the distance matrix, and use to position sets
        var distances = venn.getDistanceMatrix(areas, sets, setids);
        var positions = mds.classic(distances);

        // translate rows back to (x,y,radius) coordinates
        var circles = {};
        for (i = 0; i < sets.length; ++i) {
            var set = sets[i];
            circles[set.sets[0]] = {
                x: positions[i][0],
                y: positions[i][1],
                radius:  Math.sqrt(set.size / Math.PI)
            };
        }
        return circles;
    };

    /** Given a bunch of sets, and the desired overlaps between these sets - computes
    the distance from the actual overlaps to the desired overlaps. Note that
    this method ignores overlaps of more than 2 circles */
    venn.lossFunction = function(sets, overlaps) {
        var output = 0;

        function getCircles(indices) {
            return indices.map(function(i) { return sets[i]; });
        }

        for (var i = 0; i < overlaps.length; ++i) {
            var area = overlaps[i], overlap;
            if (area.sets.length == 1) {
                continue;
            } else if (area.sets.length == 2) {
                var left = sets[area.sets[0]],
                    right = sets[area.sets[1]];
                overlap = venn.circleOverlap(left.radius, right.radius,
                                             venn.distance(left, right));
            } else {
                overlap = venn.intersectionArea(getCircles(area.sets));
            }

            var weight = area.hasOwnProperty('weight') ? area.weight : 1.0;
            output += weight * (overlap - area.size) * (overlap - area.size);
        }

        return output;
    };

    /** Scales a solution from venn.venn or venn.greedyLayout such that it fits in
    a rectangle of width/height - with padding around the borders. also
    centers the diagram in the available space at the same time */
    venn.scaleSolution = function(solution, width, height, padding) {
        var circles = [], setids = [];
        for (var setid in solution) {
            if (solution.hasOwnProperty(setid)) {
                setids.push(setid);
                circles.push(solution[setid]);
            }
        }

        var minMax = function(d) {
            var hi = Math.max.apply(null, circles.map(
                                    function(c) { return c[d] + c.radius; } )),
                lo = Math.min.apply(null, circles.map(
                                    function(c) { return c[d] - c.radius;} ));
            return {max:hi, min:lo};
        };

        width -= 2*padding;
        height -= 2*padding;

        var xRange = minMax('x'),
            yRange = minMax('y'),
            xScaling = width  / (xRange.max - xRange.min),
            yScaling = height / (yRange.max - yRange.min),
            scaling = Math.min(yScaling, xScaling),

            // while we're at it, center the diagram too
            xOffset = (width -  (xRange.max - xRange.min) * scaling) / 2,
            yOffset = (height - (yRange.max - yRange.min) * scaling) / 2;

        var scaled = {};
        for (var i = 0; i < circles.length; ++i) {
            var circle = circles[i];
            scaled[setids[i]] = {
                radius: scaling * circle.radius,
                x: padding + xOffset + (circle.x - xRange.min) * scaling,
                y: padding + yOffset + (circle.y - yRange.min) * scaling,
            };
        }

        return scaled;
    };
})(venn);
