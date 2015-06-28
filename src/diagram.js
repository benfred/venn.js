(function(venn) {
    "use strict";
    venn.VennDiagram = function() {
        var width = 600,
            height = 350,
            padding = 15,
            duration = 1000,
            orientation = Math.PI / 2,
            normalize = true,
            fontSize = null,
            colours = d3.scale.category10(),
            layoutFunction = venn.venn;

        function chart(selection) {
            selection.each(function(data) {
                var solution = layoutFunction(data);
                if (normalize) {
                    solution = venn.normalizeSolution(solution, orientation);
                }
                var circles = venn.scaleSolution(solution, width, height, padding);
                var textCentres = computeTextCentres(circles, data);

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

        chart.normalize = function(_) {
            if (!arguments.length) return normalize;
            normalize = _;
            return chart;
        };
        chart.orientation = function(_) {
            if (!arguments.length) return orientation;
            orientation = _;
            return chart;
        };

        return chart;
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

            while (true) {
                word = words.pop();
                if (!word) break;
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

    function circleMargin(current, interior, exterior) {
        var margin = interior[0].radius - venn.distance(interior[0], current), i, m;
        for (i = 1; i < interior.length; ++i) {
            m = interior[i].radius - venn.distance(interior[i], current);
            if (m <= margin) {
                margin = m;
            }
        }

        for (i = 0; i < exterior.length; ++i) {
            m = venn.distance(exterior[i], current) - exterior[i].radius;
            if (m <= margin) {
                margin = m;
            }
        }
        return margin;
    }

    // compute the center of some circles by maximizing the margin of
    // the center point relative to the circles (interior) after subtracting
    // nearby circles (exterior)
    function computeTextCentre(interior, exterior) {
        // get an initial estimate by sampling around the interior circles
        // and taking the point with the biggest margin
        var points = [], i;
        for (i = 0; i < interior.length; ++i) {
            var c = interior[i];
            points.push({x: c.x, y: c.y});
            points.push({x: c.x + c.radius/2, y: c.y});
            points.push({x: c.x - c.radius/2, y: c.y});
            points.push({x: c.x, y: c.y + c.radius/2});
            points.push({x: c.x, y: c.y - c.radius/2});
        }
        var initial = points[0], margin = circleMargin(points[0], interior, exterior);
        for (i = 1; i < points.length; ++i) {
            var m = circleMargin(points[i], interior, exterior);
            if (m >= margin) {
                initial = points[i];
                margin = m;
            }
        }

        // maximize the margin numerically
        var solution = venn.fmin(
                    function(p) { return -1 * circleMargin({x: p[0], y: p[1]}, interior, exterior); },
                    [initial.x, initial.y],
                    {maxIterations:500, minErrorDelta:1e-10}).solution;
        var ret = {x: solution[0], y: solution[1]};

        // check solution, fallback as needed (happens if fully overlapped
        // etc)
        var valid = true;
        for (i = 0; i < interior.length; ++i) {
            if (venn.distance(ret, interior[i]) > interior[i].radius) {
                valid = false;
                break;
            }
        }

        for (i = 0; i < exterior.length; ++i) {
            if (venn.distance(ret, exterior[i]) < exterior[i].radius) {
                valid = false;
                break;
            }
        }

        if (!valid) {
            if (interior.length == 1) {
                ret = {x: interior[0].x, y: interior[0].y};
            } else {
                var areaStats = {};
                venn.intersectionArea(interior, areaStats);

                if (areaStats.arcs.length === 0) {
                    ret = {'x': 0, 'y': -1000, disjoint:true};
                } else {
                    // take average of all the points in the intersection
                    // polygon
                    ret = venn.getCenter(areaStats.arcs.map(function (a) { return a.p1; }));
                }
            }
        }

        return ret;
    }
    venn.computeTextCentre = computeTextCentre;

    function computeTextCentres(circles, areas) {
        var ret = {};
        for (var i = 0; i < areas.length; ++i) {
            var area = areas[i].sets, areaids = {};
            for (var j = 0; j < area.length; ++j) {
                areaids[area[j]] = true;
            }

            var interior = [], exterior = [];
            for (var setid in circles) {
                if (setid in areaids) {
                    interior.push(circles[setid]);
                } else {
                    exterior.push(circles[setid]);
                }
            }
            var centre = computeTextCentre(interior, exterior);
            ret[area] = centre;
            if (centre.disjoint && (areas[i].size > 0)) {
                console.log("WARNING: area " + area + " not represented on screen");
            }
        }
        return  ret;
    }
    venn.computeTextCentres = computeTextCentres;

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
})(venn);
