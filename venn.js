(function (venn) {
    "use strict";
    /** given a list of set objects, and their corresponding overlaps.
    updates the (x, y, radius) attribute on each set such that their positions
    roughly correspond to the desired overlaps */
    venn.venn = function (sets, overlaps, parameters) {
        parameters = parameters || {};
        parameters.maxIterations = parameters.maxIterations || 500;
        var lossFunction = parameters.lossFunction || venn.lossFunction;
        var initialLayout = parameters.layoutFunction || venn.greedyLayout;

        // initial layout is done greedily
        sets = initialLayout(sets, overlaps);

        // transform x/y coordinates to a vector to optimize
        var initial = new Array(2 * sets.length);
        for (var i = 0; i < sets.length; ++i) {
            initial[2 * i] = sets[i].x;
            initial[2 * i + 1] = sets[i].y;
        }

        // optimize initial layout from our loss function
        var totalFunctionCalls = 0;
        var solution = venn.fmin(
            function (values) {
                totalFunctionCalls += 1;
                var current = new Array(sets.length);
                for (var i = 0; i < sets.length; ++i) {
                    current[i] = {
                        x: values[2 * i],
                        y: values[2 * i + 1],
                        radius: sets[i].radius,
                        size: sets[i].size
                    };
                }
                return lossFunction(current, overlaps);
            },
            initial,
            parameters);

        // transform solution vector back to x/y points
        var positions = solution.solution;
        for (i = 0; i < sets.length; ++i) {
            sets[i].x = positions[2 * i];
            sets[i].y = positions[2 * i + 1];
        }

        return sets;
    };

    var SMALL = 1e-15;

    venn.circleIntegral = function (r, x) {
        var y = Math.sqrt(r * r - x * x);
        return x * y + r * r * Math.atan(x / (y + SMALL));
    };

    /** Returns the area of a circle of radius r - up to width */
    venn.circleArea = function (r, width) {
        return venn.circleIntegral(r, width - r) - venn.circleIntegral(r, -r);
    };

    /** Returns the overlap area of two circles of radius r1 and r2 - that
    have their centers separated by distance d */
    venn.circleOverlap = function (r1, r2, d) {
        // no overlap
        if (d >= r1 + r2) {
            return 0;
        }

        // completly overlapped
        if (d <= Math.abs(r1 - r2)) {
            return Math.PI * Math.min(r1, r2) * Math.min(r1, r2);
        }

        var w1 = r1 - (d * d - r2 * r2 + r1 * r1) / (2 * d),
            w2 = r2 - (d * d - r1 * r1 + r2 * r2) / (2 * d);
        return venn.circleArea(r1, w1) + venn.circleArea(r2, w2);
    };

    /** Returns the distance necessary for two circles of radius r1 + r2 to
    have the overlap area 'overlap' */
    venn.distanceFromIntersectArea = function (r1, r2, overlap) {
        if (overlap <= 0) {
            return (r1 + r2);
        }

        function loss(distance) {
            var actual;
            if (distance[0] > (r1 + r2)) {
                actual = (r1 + r2) - distance[0];
            } else {
                actual = venn.circleOverlap(r1, r2, distance[0]);
            }
            var ret = (actual - overlap) * (actual - overlap);
            return ret;
        }

        var ret = venn.fmin(loss, [Math.abs(r1 - r2)]);

        if (ret.f > 1e-3) {
            console.log("failed: " + r1 + " " + r2 + " " + overlap + " " + ret.f);
        }
        return ret.solution[0];
    };

    /// gets a matrix of euclidean distances between all sets in venn diagram
    venn.getDistanceMatrix = function (sets, overlaps) {
        // initialize an empty distance matrix between all the points
        var distances = [];
        for (var i = 0; i < sets.length; ++i) {
            distances.push([]);
            for (var j = 0; j < sets.length; ++j) {
                distances[i].push(0);
            }
        }

        // compute distances between all the points
        for (var i = 0; i < overlaps.length; ++i) {
            var current = overlaps[i];
            if (current.sets.length !== 2) {
                continue;
            }

            var left = current.sets[0],
                right = current.sets[1],
                r1 = Math.sqrt(sets[left].size / Math.PI),
                r2 = Math.sqrt(sets[right].size / Math.PI),
                distance = venn.distanceFromIntersectArea(r1, r2, current.size);
            distances[left][right] = distances[right][left] = distance;
        }
        return distances;
    };

    /** euclidean distance between two points */
    venn.distance = function (p1, p2) {
        return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                         (p1.y - p2.y) * (p1.y - p2.y));
    };

    /** Given two circles (containing a x/y/radius attributes),
    returns the intersecting points if possible.
    note: doesn't handle cases where there are infinitely many
    intersection poiints (circles are equivalent):, or only one intersection point*/
    venn.circleCircleIntersection = function (p1, p2) {
        var d = venn.distance(p1, p2),
            r1 = p1.radius,
            r2 = p2.radius;

        // if to far away, or self contained - can't be done
        if ((d >= (r1 + r2)) || (d <= Math.abs(r1 - r2))) {
            return [];
        }

        var a = (r1 * r1 - r2 * r2 + d * d) / (2 * d),
            h = Math.sqrt(r1 * r1 - a * a),
            x0 = p1.x + a * (p2.x - p1.x) / d,
            y0 = p1.y + a * (p2.y - p1.y) / d,
            rx = -(p2.y - p1.y) * (h / d),
            ry = -(p2.x - p1.x) * (h / d);

        return [{ x: x0 + rx, y: y0 - ry },
                { x: x0 - rx, y: y0 + ry }];
    };

    /** Lays out a venn diagram greedily, going from most overlapped sets to
    least overlapped, attempting to position each new set such that the
    overlapping areas to already positioned sets are basically right */
    venn.greedyLayout = function (sets, overlaps) {
        // give each set a default position + radius
        var setOverlaps = {};
        for (var i = 0; i < sets.length; ++i) {
            setOverlaps[i] = [];
            sets[i].radius = Math.sqrt(sets[i].size / Math.PI);
            sets[i].x = sets[i].y = 0;
        }

        // map each set to a list of all the other sets that overlap it
        for (i = 0; i < overlaps.length; ++i) {
            var current = overlaps[i];
            if (current.sets.length !== 2) {
                continue;
            }

            var left = current.sets[0], right = current.sets[1];
            setOverlaps[left].push({ set: right, size: current.size });
            setOverlaps[right].push({ set: left, size: current.size });
        }

        // get list of most overlapped sets
        var mostOverlapped = [];
        for (var set in setOverlaps) {
            if (setOverlaps.hasOwnProperty(set)) {
                var size = 0;
                for (i = 0; i < setOverlaps[set].length; ++i) {
                    size += setOverlaps[set][i].size;
                }

                mostOverlapped.push({ set: set, size: size });
            }
        }

        // sort by size desc
        function sortOrder(a, b) {
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
            sets[index].x = point.x;
            sets[index].y = point.y;
            positioned[index] = true;
        }

        // add most overlapped set at (0,0)
        positionSet({ x: 0, y: 0 }, mostOverlapped[0].set);

        // get distances between all points
        var distances = venn.getDistanceMatrix(sets, overlaps);

        for (i = 1; i < mostOverlapped.length; ++i) {
            var setIndex = mostOverlapped[i].set,
                set = sets[setIndex],
                overlap = setOverlaps[setIndex].filter(isPositioned);
            overlap.sort(sortOrder);

            if (overlap.length === 0) {
                throw "Need overlap information for set " + set;
            }

            var points = [];
            for (var j = 0; j < overlap.length; ++j) {
                // get appropiate distance from most overlapped already added set
                var p1 = sets[overlap[j].set],
                    d1 = distances[setIndex][overlap[j].set];

                // sample postions at 90 degrees for maximum aesheticness
                points.push({ x: p1.x + d1, y: p1.y });
                points.push({ x: p1.x - d1, y: p1.y });
                points.push({ y: p1.y + d1, x: p1.x });
                points.push({ y: p1.y - d1, x: p1.x });

                // if we have at least 2 overlaps, then figure out where the
                // set should be positioned analytically and try those too
                for (var k = j + 1; k < overlap.length; ++k) {
                    var p2 = sets[overlap[k].set],
                        d2 = distances[setIndex][overlap[k].set];

                    var extraPoints = venn.circleCircleIntersection(
                        { x: p1.x, y: p1.y, radius: d1 },
                        { x: p2.x, y: p2.y, radius: d2 });

                    for (var l = 0; l < extraPoints.length; ++l) {
                        points.push(extraPoints[l]);
                    }
                }
            }

            // we have some candidate positions for the set, examine loss
            // at each position to figure out where to put it at
            var bestLoss = 1e50, bestPoint = points[0];
            for (var j = 0; j < points.length; ++j) {
                sets[setIndex].x = points[j].x;
                sets[setIndex].y = points[j].y;
                var loss = venn.lossFunction(sets, overlaps);
                if (loss < bestLoss) {
                    bestLoss = loss;
                    bestPoint = points[j];
                }
            }

            positionSet(bestPoint, setIndex);
        }

        return sets;
    };

    /// Uses multidimensional scaling to approximate a first layout here
    venn.classicMDSLayout = function (sets, overlaps) {
        // get the distance matix
        var distances = venn.getDistanceMatrix(sets, overlaps);

        // get positions for each set
        var positions = mds.classic(distances);

        // translate back to (x,y,radius) coordinates
        for (var i = 0; i < sets.length; ++i) {
            sets[i].x = positions[i][0];
            sets[i].y = positions[i][1];
            sets[i].radius = Math.sqrt(sets[i].size / Math.PI);
        }
        return sets;
    };

    /** Given a bunch of sets, and the desired overlaps between these sets - computes
    the distance from the actual overlaps to the desired overlaps. Note that
    this method ignores overlaps of more than 2 circles */
    venn.lossFunction = function (sets, overlaps) {
        var output = 0;
        for (var i = 0; i < overlaps.length; ++i) {
            var area = overlaps[i];

            if (area.sets.length !== 2) {
                continue;
            }

            var left = sets[area.sets[0]],
                right = sets[area.sets[1]];

            var overlap = venn.circleOverlap(left.radius, right.radius,
                                             venn.distance(left, right));
            output += (overlap - area.size) * (overlap - area.size);
        }

        return output;
    };

    /** Converts an integer to a string color value */
    function intToColour(x) {
        var base = x.toString(16);
        while (base.length < 6) {
            base = "0" + base;
        }
        return "#" + base;
    }

    /** computes loss by actually laying out the sets and counting the pixels
    in a canvas. as you would expect, this is really slow. On the flipside,
    this lets us consider higher order effects.
    not only is this rather slow - it also may have issues. I'm not convinced
    about anti-aliasing is truly disabled here/1*/
    venn.canvasLossFunction = function (sets, overlaps) {
        var canvas = document.createElement("canvas");
        canvas.width = canvas.height = 200;

        var context = canvas.getContext('2d');
        context.globalCompositeOperation = 'lighter';
        context.webkitImageSmoothingEnabled = false;

        var scaled = venn.scaleSolution(sets, canvas.width, canvas.height, 0),
            scaling = scaled.scaling;

        for (var i = 0; i < scaled.length; ++i) {
            var set = scaled[i];
            context.beginPath();
            context.arc(set.x, set.y, set.radius, 0, 2 * Math.PI, false);
            context.fillStyle = intToColour(1 << i);
            context.closePath();
            context.fill();
        }

        var histogram = {},
            data = context.getImageData(0, 0, canvas.width, canvas.height).data;
        for (i = 0; i < data.length; i += 4) {
            var colour = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
            if (colour in histogram) {
                histogram[colour] += 1;
            } else {
                histogram[colour] = 1;
            }
        }

        var out = 0;
        for (i = 0; i < overlaps.length; ++i) {
            var key = 0;
            for (var j = 0; j < overlaps[i].sets.length; j++) {
                key |= (1 << overlaps[i].sets[j]);
            }

            var desiredPixels = overlaps[i].size * scaling * scaling;
            var actualPixels = 0;
            for (var other in histogram) {
                if (histogram.hasOwnProperty(other)) {
                    var otherKey = parseInt(other, 10);
                    if ((otherKey & key) === key) {
                        actualPixels += histogram[other];
                    }
                }
            }

            var delta = (desiredPixels - actualPixels);
            out += delta * delta;
        }

        return out;
    };

    /** Scales a solution from venn.venn or venn.greedyLayout such that it fits in
    a rectangle of width/height - with padding around the borders. */
    venn.scaleSolution = function (solution, width, height, padding) {
        var minMax = function (d) {
            var hi = Math.max.apply(null, solution.map(
                                    function (c) { return c[d] + c.radius; })),
                lo = Math.min.apply(null, solution.map(
                                    function (c) { return c[d] - c.radius; }));
            return { max: hi, min: lo };
        };

        width -= 2 * padding;
        height -= 2 * padding;

        var xRange = minMax('x'),
            yRange = minMax('y'),
            xScaling = width / (xRange.max - xRange.min),
            yScaling = height / (yRange.max - yRange.min),
            scaling = Math.min(yScaling, xScaling);

        for (var i = 0; i < solution.length; ++i) {
            var set = solution[i];
            set.radius = scaling * set.radius;
            set.x = padding + (set.x - xRange.min) * scaling;
            set.y = padding + (set.y - yRange.min) * scaling;
        }
        solution.scaling = scaling;

        return solution;
    };

    function weightedSum(a, b) {
        var ret = new Array(a[1].length || 0);
        for (var j = 0; j < ret.length; ++j) {
            ret[j] = a[0] * a[1][j] + b[0] * b[1][j];
        }
        return ret;
    }

    /** minimizes a function using the downhill simplex method */
    venn.fmin = function (f, x0, parameters) {
        parameters = parameters || {};

        var maxIterations = parameters.maxIterations || x0.length * 200,
            nonZeroDelta = parameters.nonZeroDelta || 1.1,
            zeroDelta = parameters.zeroDelta || 0.001,
            minErrorDelta = parameters.minErrorDelta || 1e-5,
            rho = parameters.rho || 1,
            chi = parameters.chi || 2,
            psi = parameters.psi || -0.5,
            sigma = parameters.sigma || 0.5,
            callback = parameters.callback;

        // initialize simplex.
        var N = x0.length,
            simplex = new Array(N + 1);
        simplex[0] = x0;
        simplex[0].fx = f(x0);
        for (var i = 0; i < N; ++i) {
            var point = x0.slice();
            point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
            simplex[i + 1] = point;
            simplex[i + 1].fx = f(point);
        }

        var sortOrder = function (a, b) { return a.fx - b.fx; };

        for (var iteration = 0; iteration < maxIterations; ++iteration) {
            simplex.sort(sortOrder);
            if (callback) {
                callback(simplex);
            }

            if (Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) {
                break;
            }

            // compute the centroid of all but the worst point in the simplex
            var centroid = new Array(N);
            for (i = 0; i < N; ++i) {
                centroid[i] = 0;
                for (var j = 0; j < N; ++j) {
                    centroid[i] += simplex[j][i];
                }
                centroid[i] /= N;
            }

            // reflect the worst point past the centroid  and compute loss at reflected
            // point
            var worst = simplex[N];
            var reflected = weightedSum([1 + rho, centroid], [-rho, worst]);
            reflected.fx = f(reflected);

            var replacement = reflected;

            // if the reflected point is the best seen, then possibly expand
            if (reflected.fx <= simplex[0].fx) {
                var expanded = weightedSum([1 + chi, centroid], [-chi, worst]);
                expanded.fx = f(expanded);
                if (expanded.fx < reflected.fx) {
                    replacement = expanded;
                }
            }

                // if the reflected point is worse than the second worst, we need to
                // contract
            else if (reflected.fx >= simplex[N - 1].fx) {
                var shouldReduce = false;
                var contracted;

                if (reflected.fx <= worst.fx) {
                    // do an inside contraction
                    contracted = weightedSum([1 + psi, centroid], [-psi, worst]);
                    contracted.fx = f(contracted);
                    if (contracted.fx < worst.fx) {
                        replacement = contracted;
                    } else {
                        shouldReduce = true;
                    }
                } else {
                    // do an outside contraction
                    contracted = weightedSum([1 - psi * rho, centroid], [psi * rho, worst]);
                    contracted.fx = f(contracted);
                    if (contracted.fx <= reflected.fx) {
                        replacement = contracted;
                    } else {
                        shouldReduce = true;
                    }
                }

                if (shouldReduce) {
                    // do reduction. doesn't actually happen that often
                    for (i = 1; i < simplex.length; ++i) {
                        simplex[i] = weightedSum([1 - sigma, simplex[0]],
                                                 [sigma - 1, simplex[i]]);
                        simplex[i].fx = f(simplex[i]);
                    }
                }
            }

            simplex[N] = replacement;
        }

        simplex.sort(sortOrder);
        return {
            f: simplex[0].fx,
            solution: simplex[0]
        };
    };

    venn.drawD3Diagram = function (element, dataset, width, height, padding) {
        var pi = Math.PI;
        padding = padding || 6;

        dataset = venn.scaleSolution(dataset, width, height, padding);

        var svg = element.append("svg")
                .attr("width", width)
                .attr("height", height);

        var colours = d3.scale.category10();

        var nodes = svg.selectAll("circle")
                         .data(dataset)
                         .enter()
                         .append("g");

        nodes.append("circle")
               .attr("r", function (d) { return d.radius; })
               .style("fill-opacity", 0.2)
               .style("stroke-opacity", 0.8)
               .style("stroke-width", 3)
               .attr("cx", function (d) { return d.x; })
               .attr("cy", function (d) { return d.y; })
               .style("stroke", function (d, i) { return colours(i); })
               .style("fill", function (d, i) { return colours(i); });

        var arcs = d3.svg.arc()
                .innerRadius(function (d) { return d.radius }) // make really small 
                .outerRadius(function (d) { return d.radius })
                .startAngle(-pi)
                .endAngle(pi);

        nodes.append("path")
                .attr("d", arcs)
                .attr("id", function(d, i) { return "path" + i; })
                .attr("transform", function (d) { return "translate(" + d.x + ", " + d.y + ")"; })
                .style("stroke-opacity", 0.5)
                .style("fill", function (d, i) { return colours(i); });


        // formats number with commas
        function FormatNumber(yourNumber) {
            //Seperates the components of the number
            var n = yourNumber.toString().split(".");
            //Comma-fies the first part
            n[0] = n[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            //Combines the two sections
            return n.join(".");
        }
        
        // display the size of the circle in the middle
        nodes.append("text")
               .attr("x", function (d) { return d.x; })
               .attr("y", function (d) { return d.y; })
               .attr("text-anchor", "middle")
               .style("stroke", function (d, i) { return colours(i); })
               .style("fill", function (d, i) { return colours(i); })
               .text(function (d) { return FormatNumber(d.size); });
        
        // display the title of the circle along the path
        nodes.append("text")
                .attr("x", function (d) { return d.x })
                .attr("dy", 15)
                .append("textPath")
                .style("stroke", function (d, i) { return colours(i); })
                .style("fill", function (d, i) { return colours(i); })
                .attr("xlink:href", function (d, i) { return "#path" + i })
                .attr("letter-spacing", "0.15em")
                .text(function (d) { return d.label });

        var colors = new Array();
        for (var i = 0; i < dataset.length; i++) {
            colors.push(colours(i));
        }

        // return array of colors used to be later used in the 
        // legend if needed
        return colors;
    };

    venn.updateD3Diagram = function (element, dataset) {
        var svg = element.select("svg"),
            width = parseInt(svg.attr('width'), 10),
            height = parseInt(svg.attr('height'), 10);

        dataset = venn.scaleSolution(dataset, width, height, 6);
        element.selectAll("circle")
            .data(dataset)
            .transition()
            .duration(400)
            .attr("cx", function (d) { return d.x; })
            .attr("cy", function (d) { return d.y; })
            .attr("r", function (d) { return d.radius; });

        element.selectAll("text")
            .data(dataset)
            .transition()
            .duration(400)
            .attr("x", function (d) { return d.x; })
            .attr("y", function (d) { return d.y; });
    };
}(window.venn = window.venn || {}));
