 (function(venn) {
    venn.VennDiagram = function() {
        var width = 600,
            height = 350,
            padding = 15,
            duration = 1000,
            fontSize = null,
            colours = d3.scale.category10(),
            layoutFunction = venn.venn;

        function chart(selection) {
            selection.each(function(data) {
                // calculate circle position, scale to fit
                var circles = venn.scaleSolution(layoutFunction(data), width, height, padding);
                var textCentres = computeTextCentres(circles, data, width, height);

                // draw out a svg
                var svg = d3.select(this).selectAll("svg").data([circles]);
                svg.enter().append("svg");

                svg.attr("width", width)
                   .attr("height", height);

                // to properly transition intersection areas, we need the
                // previous circles locations. load from elements
                var previous = {}, hasPrevious = false;
                svg.selectAll("g").each(function (d) {
                    var path = d3.select(this).select("path").attr("d");
                    if ((d.sets.length == 1) && path) {
                        hasPrevious = true;
                        previous[d.sets[0]] = venn.circleFromPath(path);
                    }
                });

                // interpolate intersection area paths between previous and
                // current paths
                var pathTween = function(d) {
                    return function(t) {
                        var c = d.sets.map(function (set) {
                            var start = previous[set], end = circles[set];
                            if (!start) {
                                start = {x : width/2, y : height/2, radius : 1};
                            }
                            if (!end) {
                                end = {x : width/2, y : height/2, radius : 1};
                            }
                            return {'x' : start.x * (1 - t) + end.x * t,
                                    'y' : start.y * (1 - t) + end.y * t,
                                    'radius' : start.radius * (1 - t) + end.radius * t};

                        });
                        return venn.intersectionAreaPath(c);
                    };
                };

                // update data, joining on the set ids
                var nodes = svg.selectAll("g")
                    .data(data, function(d) { return d.sets; });

                // create new nodes
                var enter = nodes.enter()
                    .append('g')
                    .attr("class", function(d) {
                        return "venn-area venn-" +
                            (d.sets.length == 1 ? "circle" : "intersection") +
                            (" venn-sets-" + d.sets.join("_"));
                    });

                enter.append("path")
                    .style("fill-opacity", "0")
                    .filter(function (d) { return d.sets.length == 1; } )
                    .style("fill", function(d) { return colours(label(d)); })
                    .style("fill-opacity", ".25");

                var enterText = enter.append("text")
                    .style("fill", function(d) { return d.sets.length == 1 ? colours(label(d)) : "#444"; })
                    .text(function (d) { return label(d); } )
                    .attr("text-anchor", "middle")
                    .attr("dy", ".35em")
                    .attr("x", width/2)
                    .attr("y", height/2);

                // update existing
                var update = nodes.transition("venn").duration(hasPrevious ? duration : 0);
                update.select("path")
                    .attrTween("d", pathTween);

                var updateText = update.select("text")
                    .text(function (d) { return label(d); } )
                    .each("end", venn.wrapText(circles, label))
                    .attr("x", function(d) {
                        return Math.floor(textCentres[d.sets].x);
                    })
                    .attr("y", function(d) {
                        return Math.floor(textCentres[d.sets].y);
                    });

                // if we've been passed a fontSize explicitly, use it to
                // transition
                if (fontSize !== null) {
                    enterText.style("font-size", "0px");
                    updateText.style("font-size", fontSize);
                }

                // remove old
                var remove = nodes.exit().transition('venn').duration(duration).remove();
                remove.select("path")
                    .attrTween("d", pathTween);

                remove.select("text")
                    .text(function (d) { return label(d); } )
                    .attr("x", width/2)
                    .attr("y", height/2)
                    .style("font-size", "0px");
            });
        }

        function label(d) {
            if (d.label) {
                return d.label;
            }
            if (d.sets.length == 1) {
                return '' + d.sets[0];
            }
        }

        chart.width = function(_) {
            if (!arguments.length) return width;
            width = _;
            return chart;
        };

        chart.height = function(_) {
            if (!arguments.length) return height;
            height = _;
            return chart;
        };

        chart.padding = function(_) {
            if (!arguments.length) return padding;
            padding = _;
            return chart;
        };

        chart.colours = function(_) {
            if (!arguments.length) return colours;
            colours = _;
            return chart;
        };

        chart.fontSize = function(_) {
            if (!arguments.length) return fontSize;
            fontSize = _;
            return chart;
        };

        chart.duration = function(_) {
            if (!arguments.length) return duration;
            duration = _;
            return chart;
        };

        chart.layoutFunction = function(_) {
            if (!arguments.length) return layoutFunction;
            layoutFunction = _;
            return chart;
        };

        return chart;
    };

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

    // sometimes text doesn't fit inside the circle, if thats the case lets wrap
    // the text here such that it fits
    // todo: looks like this might be merged into d3 (
    // https://github.com/mbostock/d3/issues/1642),
    // also worth checking out is
    // http://engineering.findthebest.com/wrapping-axis-labels-in-d3-js/
    // this seems to be one of those things that should be easy but isn't
    venn.wrapText = function(circles, labeller) {
        return function() {
            var text = d3.select(this),
                data = text.datum(),
                width = circles[data.sets[0]].radius || 50,
                label = labeller(data) || '';

                var words = label.split(/\s+/).reverse(),
                maxLines = 3,
                minChars = (label.length + words.length) / maxLines,
                word = words.pop(),
                line = [word],
                joined,
                lineNumber = 0,
                lineHeight = 1.1, // ems
                tspan = text.text(null).append("tspan").text(word);

            while (word = words.pop()) {
                line.push(word);
                joined = line.join(" ");
                tspan.text(joined);
                if (joined.length > minChars && tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").text(word);
                    lineNumber++;
                }
            }

            var initial = 0.35 - lineNumber * lineHeight / 2,
                x = text.attr("x"),
                y = text.attr("y");

            text.selectAll("tspan")
                .attr("x", x)
                .attr("y", y)
                .attr("dy", function(d, i) {
                     return (initial + i * lineHeight) + "em";
                });
        };
    };

    /** finds the zeros of a function, given two starting points (which must
     * have opposite signs */
    venn.bisect = function(f, a, b, parameters) {
        parameters = parameters || {};
        var maxIterations = parameters.maxIterations || 100,
            tolerance = parameters.tolerance || 1e-10,
            fA = f(a),
            fB = f(b),
            delta = b - a;

        if (fA * fB > 0) {
            throw "Initial bisect points must have opposite signs";
        }

        if (fA === 0) return a;
        if (fB === 0) return b;

        for (var i = 0; i < maxIterations; ++i) {
            delta /= 2;
            var mid = a + delta,
                fMid = f(mid);

            if (fMid * fA >= 0) {
                a = mid;
            }

            if ((Math.abs(delta) < tolerance) || (fMid === 0)) {
                return mid;
            }
        }
        return a + delta;
    };

    function weightedSum(ret, w1, v1, w2, v2) {
        for (var j = 0; j < ret.length; ++j) {
            ret[j] = w1 * v1[j] + w2 * v2[j];
        }
    }


    /** minimizes a function using the downhill simplex method */
    venn.fmin = function(f, x0, parameters) {
        parameters = parameters || {};

        var maxIterations = parameters.maxIterations || x0.length * 200,
            nonZeroDelta = parameters.nonZeroDelta || 1.1,
            zeroDelta = parameters.zeroDelta || 0.001,
            minErrorDelta = parameters.minErrorDelta || 1e-6,
            rho = parameters.rho || 1,
            chi = parameters.chi || 2,
            psi = parameters.psi || -0.5,
            sigma = parameters.sigma || 0.5,
            callback = parameters.callback,
            temp;

        // initialize simplex.
        var N = x0.length,
            simplex = new Array(N + 1);
        simplex[0] = x0;
        simplex[0].fx = f(x0);
        for (var i = 0; i < N; ++i) {
            var point = x0.slice();
            point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
            simplex[i+1] = point;
            simplex[i+1].fx = f(point);
        }

        var sortOrder = function(a, b) { return a.fx - b.fx; };

        var centroid = x0.slice(),
            reflected = x0.slice(),
            contracted = x0.slice(),
            expanded = x0.slice();

        for (var iteration = 0; iteration < maxIterations; ++iteration) {
            simplex.sort(sortOrder);
            if (callback) {
                callback(simplex);
            }

            if (Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) {
                break;
            }

            // compute the centroid of all but the worst point in the simplex
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
            weightedSum(reflected, 1+rho, centroid, -rho, worst);
            reflected.fx = f(reflected);

            // if the reflected point is the best seen, then possibly expand
            if (reflected.fx <= simplex[0].fx) {
                weightedSum(expanded, 1+chi, centroid, -chi, worst);
                expanded.fx = f(expanded);
                if (expanded.fx < reflected.fx) {
                    temp = simplex[N];
                    simplex[N] = expanded;
                    expanded = temp;
                }  else {
                    temp = simplex[N];
                    simplex[N] = reflected;
                    reflected = temp;
                }
            }

            // if the reflected point is worse than the second worst, we need to
            // contract
            else if (reflected.fx >= simplex[N-1].fx) {
                var shouldReduce = false;

                if (reflected.fx <= worst.fx) {
                    // do an inside contraction
                    weightedSum(contracted, 1+psi, centroid, -psi, worst);
                    contracted.fx = f(contracted);
                    if (contracted.fx < worst.fx) {
                        temp = simplex[N];
                        simplex[N] = contracted;
                        contracted = temp;
                    } else {
                        shouldReduce = true;
                    }
                } else {
                    // do an outside contraction
                    weightedSum(contracted, 1-psi * rho, centroid, psi*rho, worst);
                    contracted.fx = f(contracted);
                    if (contracted.fx <= reflected.fx) {
                        temp = simplex[N];
                        simplex[N] = contracted;
                        contracted = temp;
                    } else {
                        shouldReduce = true;
                    }
                }

                if (shouldReduce) {
                    // do reduction. doesn't actually happen that often
                    for (i = 1; i < simplex.length; ++i) {
                        weightedSum(simplex[i], 1 - sigma, simplex[0], sigma - 1, simplex[i]);
                        simplex[i].fx = f(simplex[i]);
                    }
                }
            } else {
                temp = simplex[N];
                simplex[N] = reflected;
                reflected = temp;
            }

        }

        simplex.sort(sortOrder);
        return {f : simplex[0].fx,
                solution : simplex[0]};
    };

    venn.circlePath = function(x, y, r) {
        var ret = [];
        ret.push("\nM", x, y);
        ret.push("\nm", -r, 0);
        ret.push("\na", r, r, 0, 1, 0, r *2, 0);
        ret.push("\na", r, r, 0, 1, 0,-r *2, 0);
        return ret.join(" ");
    };

    // inverse of the circlePath function, returns a circle object from an svg path
    venn.circleFromPath = function(path) {
        var tokens = path.split(' ');
        return {'x' : parseFloat(tokens[1]),
                'y' : parseFloat(tokens[2]),
                'radius' : -parseFloat(tokens[4])
                };
    };

    /** returns a svg path of the intersection area of a bunch of circles */
    venn.intersectionAreaPath = function(circles) {
        var stats = {};
        venn.intersectionArea(circles, stats);
        var arcs = stats.arcs;

        if (arcs.length === 0) {
            return "M 0 0";

        } else if (arcs.length == 1) {
            var circle = arcs[0].circle;
            return venn.circlePath(circle.x, circle.y, circle.radius);

        } else {
            // draw path around arcs
            var ret = ["\nM", arcs[0].p2.x, arcs[0].p2.y];
            for (var i = 0; i < arcs.length; ++i) {
                var arc = arcs[i], r = arc.circle.radius, wide = arc.width > r;
                ret.push("\nA", r, r, 0, wide ? 1 : 0, 1,
                         arc.p1.x, arc.p1.y);
            }
            return ret.join(" ");
        }
    };

    function computeTextCentres(circles, areas, width, height) {
        // basically just finding the center point of each region by sampling
        // points in a grid.
        var points = {};

        var samples = 32;
        for (var i = 0; i < samples; ++i) {
            var x = i * width / samples;
            for (var j = 0; j < samples; ++j) {
                var y = j * height / samples;
                var point = {'x' : x, 'y' : y};

                var contained = [];
                for (var k in circles) {
                    if (venn.distance(point, circles[k]) <= circles[k].radius) {
                        contained.push(k);
                    }
                }
                if (!contained.length) {
                    continue;
                }

                if (!(contained in points)) {
                    points[contained] = [];
                }
                points[contained].push(point);
            }
        }

        function getCircles(area) {
            return area.map(function(a) { return circles[a]; });
        }

        var ret = {};
        for (i = 0; i < areas.length; ++i) {
            var area = areas[i].sets;
            if (points.hasOwnProperty(area)) {
                ret[area] = venn.getCenter(points[area]);
                // todo: if the point isn't in the area (could not even be a
                // circle), find the nearest point that is in the region,
                // extend past that
            } else if (area.length == 1) {
                // small/missing intersection area. default to circle centre
                var circle = circles[area[0]];
                ret[area] = {'x': circle.x, 'y': circle.y };

            } else {
                var areaStats = {};
                venn.intersectionArea(getCircles(area), areaStats);

                if (areaStats.arcs.length === 0) {
                    if (areas[i].size > 0) {
                        console.log("WARNING: area " + JSON.stringify(area) + " isn't represented on diagram" );
                    }
                    ret[area] = {'x' : 0, 'y' : -1000};
                } else {
                    // take average of all the points in the outer perimiter
                    ret[area] = venn.getCenter(areaStats.arcs.map(function (a) { return a.p1; }));
                }
            }
        }
        return ret;
    }

    // sorts all areas in the venn diagram, so that
    // a particular area is on top (relativeTo) - and
    // all other areas are so that the smallest areas are on top
    venn.sortAreas = function(div, relativeTo) {
        // need to sort div's so that Z order is correct
        div.selectAll("g").sort(function (a, b) {
            // highest order set intersections first
            if (a.sets.length != b.sets.length) {
                return a.sets.length - b.sets.length;
            }

            // current element is highest inside its order
            if ((a == relativeTo) || (b == relativeTo)) {
                return (a == relativeTo) ? 1 : -1;
            }

            // finally by size
            return b.size - a.size;
        });
    }


    var SMALL = 1e-10;

    /** Returns the intersection area of a bunch of circles (where each circle
     is an object having an x,y and radius property) */
    venn.intersectionArea = function(circles, stats) {
        // get all the intersection points of the circles
        var intersectionPoints = getIntersectionPoints(circles);

        // filter out points that aren't included in all the circles
        var innerPoints = intersectionPoints.filter(function (p) {
            return venn.containedInCircles(p, circles);
        });

        var arcArea = 0, polygonArea = 0, arcs = [], i;

        // if we have intersection points that are within all the circles,
        // then figure out the area contained by them
        if (innerPoints.length > 1) {
            // sort the points by angle from the center of the polygon, which lets
            // us just iterate over points to get the edges
            var center = venn.getCenter(innerPoints);
            for (i = 0; i < innerPoints.length; ++i ) {
                var p = innerPoints[i];
                p.angle = Math.atan2(p.x - center.x, p.y - center.y);
            }
            innerPoints.sort(function(a,b) { return b.angle - a.angle;});

            // iterate over all points, get arc between the points
            // and update the areas
            var p2 = innerPoints[innerPoints.length - 1];
            for (i = 0; i < innerPoints.length; ++i) {
                var p1 = innerPoints[i];

                // polygon area updates easily ...
                polygonArea += (p2.x + p1.x) * (p1.y - p2.y);

                // updating the arc area is a little more involved
                var midPoint = {x : (p1.x + p2.x) / 2,
                                y : (p1.y + p2.y) / 2},
                    arc = null;

                for (var j = 0; j < p1.parentIndex.length; ++j) {
                    if (p2.parentIndex.indexOf(p1.parentIndex[j]) > -1) {
                        // figure out the angle halfway between the two points
                        // on the current circle
                        var circle = circles[p1.parentIndex[j]],
                            a1 = Math.atan2(p1.x - circle.x, p1.y - circle.y),
                            a2 = Math.atan2(p2.x - circle.x, p2.y - circle.y);

                        var angleDiff = (a2 - a1);
                        if (angleDiff < 0) {
                            angleDiff += 2*Math.PI;
                        }

                        // and use that angle to figure out the width of the
                        // arc
                        var a = a2 - angleDiff/2,
                            width = venn.distance(midPoint, {
                                x : circle.x + circle.radius * Math.sin(a),
                                y : circle.y + circle.radius * Math.cos(a)
                            });

                        // pick the circle whose arc has the smallest width
                        if ((arc === null) || (arc.width > width)) {
                            arc = { circle : circle,
                                    width : width,
                                    p1 : p1,
                                    p2 : p2};
                        }
                    }
                }
                arcs.push(arc);
                arcArea += venn.circleArea(arc.circle.radius, arc.width);
                p2 = p1;
            }
        } else {
            // no intersection points, is either disjoint - or is completely
            // overlapped. figure out which by examining the smallest circle
            var smallest = circles[0];
            for (i = 1; i < circles.length; ++i) {
                if (circles[i].radius < smallest.radius) {
                    smallest = circles[i];
                }
            }

            // make sure the smallest circle is completely contained in all
            // the other circles
            var disjoint = false;
            for (i = 0; i < circles.length; ++i) {
                if (venn.distance(circles[i], smallest) > Math.abs(smallest.radius - circles[i].radius)) {
                    disjoint = true;
                    break;
                }
            }

            if (disjoint) {
                arcArea = polygonArea = 0;

            } else {
                arcArea = smallest.radius * smallest.radius * Math.PI;
                arcs.push({circle : smallest,
                           p1: { x: smallest.x,        y : smallest.y + smallest.radius},
                           p2: { x: smallest.x - SMALL, y : smallest.y + smallest.radius},
                           width : smallest.radius * 2 });
            }
        }

        polygonArea /= 2;
        if (stats) {
            stats.area = arcArea + polygonArea;
            stats.arcArea = arcArea;
            stats.polygonArea = polygonArea;
            stats.arcs = arcs;
            stats.innerPoints = innerPoints;
            stats.intersectionPoints = intersectionPoints;
        }

        return arcArea + polygonArea;
    };

    /** returns whether a point is contained by all of a list of circles */
    venn.containedInCircles = function(point, circles) {
        for (var i = 0; i < circles.length; ++i) {
            if (venn.distance(point, circles[i]) > circles[i].radius + SMALL) {
                return false;
            }
        }
        return true;
    };

    /** Gets all intersection points between a bunch of circles */
    function getIntersectionPoints(circles) {
        var ret = [];
        for (var i = 0; i < circles.length; ++i) {
            for (var j = i + 1; j < circles.length; ++j) {
                var intersect = venn.circleCircleIntersection(circles[i],
                                                              circles[j]);
                for (var k = 0; k < intersect.length; ++k) {
                    var p = intersect[k];
                    p.parentIndex = [i,j];
                    ret.push(p);
                }
            }
        }
        return ret;
    }

    venn.circleIntegral = function(r, x) {
        var y = Math.sqrt(r * r - x * x);
        return x * y + r * r * Math.atan2(x, y);
    };

    /** Returns the area of a circle of radius r - up to width */
    venn.circleArea = function(r, width) {
        return venn.circleIntegral(r, width - r) - venn.circleIntegral(r, -r);
    };


    /** euclidean distance between two points */
    venn.distance = function(p1, p2) {
        return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                         (p1.y - p2.y) * (p1.y - p2.y));
    };


    /** Returns the overlap area of two circles of radius r1 and r2 - that
    have their centers separated by distance d. Simpler faster
    circle intersection for only two circles */
    venn.circleOverlap = function(r1, r2, d) {
        // no overlap
        if (d >= r1 + r2) {
            return 0;
        }

        // completely overlapped
        if (d <= Math.abs(r1 - r2)) {
            return Math.PI * Math.min(r1, r2) * Math.min(r1, r2);
        }

        var w1 = r1 - (d * d - r2 * r2 + r1 * r1) / (2 * d),
            w2 = r2 - (d * d - r1 * r1 + r2 * r2) / (2 * d);
        return venn.circleArea(r1, w1) + venn.circleArea(r2, w2);
    };


    /** Given two circles (containing a x/y/radius attributes),
    returns the intersecting points if possible.
    note: doesn't handle cases where there are infinitely many
    intersection points (circles are equivalent):, or only one intersection point*/
    venn.circleCircleIntersection = function(p1, p2) {
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

        return [{ x: x0 + rx, y : y0 - ry },
                { x: x0 - rx, y : y0 + ry }];
    };

    /** Returns the center of a bunch of points */
    venn.getCenter = function(points) {
        var center = { x: 0, y: 0};
        for (var i =0; i < points.length; ++i ) {
            center.x += points[i].x;
            center.y += points[i].y;
        }
        center.x /= points.length;
        center.y /= points.length;
        return center;
    };
}(window.venn = window.venn || {}));
