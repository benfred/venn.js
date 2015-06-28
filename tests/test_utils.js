function randomPoint(rect) {
    return {x: rect.x + Math.random() * rect.width,
            y: rect.y + Math.random() * rect.height};
}

function generateRandomCircles(count, minRadius, maxRadius) {
    var range = {x : 0,  y : 0, height : 1, width : 1},
        output = [];

    for (var i = 0; i < count; ++i ) {
        var p = randomPoint(range);
        p.radius = Math.random() * (maxRadius - minRadius) + minRadius;
        output.push(p);
    }

    return output;
}

function getCombinations(input) {
    var output = [];
    function inner(i, current) {
        for (var j = i + 1; j < input.length; ++j) {
            var next = current.slice();
            next.push(input[j]);
            output.push(next);
            inner(j, next);
        }
    }
    inner(-1, []);
    return output;
}

function getAllIntersections(circles) {
    var setids = circles.map(function (x, i) { return i; });
    return getCombinations(setids).map(function (ids) {
        return {sets: ids,
                size: venn.intersectionArea(ids.map(function(i) {
                    return circles[i]; }))
            };
    });
}

function getIntersections(circles) {
    return getAllIntersections(circles)
        .filter(function(a) { return (a.sets.length <= 2); });
}

function getErrors(areas, current) {
    var intersections = areas.filter(function (a) {
        return a.sets.length ==2;
    });

    var total = 0, i;
    for (i = 0; i < areas.length; ++i) {
        total += areas[i].size;
    }

    var ret = [];
    for (i = 0; i < areas.length; ++i) {
        var area = areas[i];
        if (area.sets.length == 2) {
            var left = current[area.sets[0]], right = current[area.sets[1]],
                overlap = venn.circleOverlap(left.radius, right.radius,
                                            venn.distance(left, right)),
                delta = Math.abs(overlap - area.size);

            if (isNaN(delta)) {
                ret.push(area);
            }

            else if ((delta / (1e-10 + Math.max(area.size, overlap)) >= 0.1) &&
                (delta > total / (50000))) {
                ret.push(area);
            }
        }
    }
    return ret;
}

function hilightErrors(div, areas, current, duration) {
    // calculate success, and highlight failures in red
    var failedAreas = getErrors(areas, current);

    div.selectAll(".venn-intersection path")
        .style("stroke-width", 0).style("fill-opacity", 0);

    for (var i = 0; i < failedAreas.length; ++i) {
        var area = failedAreas[i];
        div.selectAll(".venn-sets-" + area.sets[0] + "_" + area.sets[1] + " path")
            .style("stroke", "red")
            .style("stroke-opacity", 1)
            .style("stroke-width", 2)
            .style("fill", "red")
            .style("fill-opacity", 0.2);
    }
    return failedAreas.length === 0;
}

// converts our areas into disjoint areas for use by other libraries
// like matplotlib_venn, venneuler etc
function makeDisjoint(areas) {
    var lookup = {};
    var ret = [], i;
    for (i = 0; i < areas.length; ++i) {
        var clone = {sets: areas[i].sets, size:areas[i].size};
        lookup[areas[i].sets] = clone;
        ret.push(clone);
    }

    ret.sort(function(a, b) {
        return b.sets.length - a.sets.length;
    });

    for (i = 0; i < ret.length; ++i) {
        var area = ret[i];
        var supersets = getCombinations(area.sets);
        for (var j = 0; j < supersets.length; ++j) {
            if (supersets[j].length < area.sets.length) {
                lookup[supersets[j]].size -= area.size;
            }
        }
    }
    return ret;
}

var MIN_CIRCLES = 2, MAX_CIRCLES = 8;
function getCircleCategories() {
    var ret = [];
    for (var i = MIN_CIRCLES; i <= MAX_CIRCLES; ++i) {
        ret.push(i);
    }
    return ret;
}

// tracks status of perfomance of each algorithm
function VennAlgorithm(name, method) {
    this.name = name;
    this.method = method;
    this.errors = {};
    this.counts = {};
    this.times = {};

    for (var i = MIN_CIRCLES; i <= MAX_CIRCLES; ++i) {
        this.errors[i] = 0;
        this.counts[i] = 0;
        this.times[i] = 0;
    }
}

// test out the algorithm with the current circles/intersection areas
VennAlgorithm.prototype.update = function(intersections, computed, elapsed) {
    var count = intersections.filter(function(a) { return a.sets.length == 1; }).length;
    if (getErrors(intersections, computed).length > 0)  {
        this.errors[count] += 1;
    }
    this.counts[count] += 1;
    this.times[count] += elapsed;
};

VennAlgorithm.prototype.reset = function() {
    for (var i = MIN_CIRCLES; i <= MAX_CIRCLES; ++i) {
        this.errors[i] = 0;
        this.counts[i] = 0;
        this.times[i] = 0;
    }
};

// random layout: poor choice
function randomLayout(areas) {
    var range = {x : 0,  y : 0, height : 1, width : 1};
    var ret = [];
    for (var i = 0; i < areas.length; ++i) {
        var area = areas[i];
        if (area.sets.length == 1) {
            ret.push({x : Math.random(),
                      y : Math.random(),
                      radius: Math.sqrt(area.size / Math.PI)});
        }
    }
    return ret;
}

function bestInitial(areas) {
    var g = venn.greedyLayout(areas);
    var w  = constrainedMDSLayout(areas, {restarts:5});
    if (venn.lossFunction(w, areas) < venn.lossFunction(g, areas)) {
        g = w;
    }
    return g;
}

function createPerformanceChart(selector,
                                algorithms,
                                precalcPerf,
                                precalcSpeed) {
    function getPerformanceData() {
        var columns = [];
        for (var i = 0; i < algorithms.length; ++i) {
            var algo = algorithms[i];
            var column = [algo.name];
            for (var j = MIN_CIRCLES; j <= MAX_CIRCLES; ++j) {

                var c = algo.counts[j] || 1;
                column.push((c - algo.errors[j] ) / c);
            }
            columns.push(column);
        }

        return { 'columns' : columns };
    }
    this.getPerformanceData = getPerformanceData;

    function getSpeedData() {
        var columns = [];
        for (var i = 0; i < algorithms.length; ++i) {
            var algo = algorithms[i];
            var column = [algo.name];
            for (var j = MIN_CIRCLES; j <= MAX_CIRCLES; ++j) {
                var c = algo.counts[j] || 1;
                column.push(algo.times[j] / c);
            }
            columns.push(column);
        }

        return { 'columns' : columns };
    }

    this.getSpeedData = getSpeedData;

    var performanceChart = c3.generate({
        bindto: selector + " #perfchart",
        axis: {
            x: {
                label: {text: 'Number of Sets', position: 'outer-center'},
                type: 'category',
                categories: getCircleCategories(),
            },
            y: {label: {text: 'Success Rate', position: 'outer-middle'},
                min:0.05,
                max:1,
                tick : {format: function(d) { return (100*d).toFixed(0) +
                "%";}}
            },
        },
        tooltip : {
            format : {
                title: function(x) { return (x + MIN_CIRCLES) + " Sets";} ,
                value: function(x) { return "" + (100*x).toFixed(2) + "%"; },
            },
        },
        data:precalcPerf || getPerformanceData() });


    var speedChart = c3.generate({
        bindto: selector + " #speedchart",
        axis: {
            x: {
                label: {text: 'Number of Sets', position: 'outer-center'},
                type: 'category',
                categories: getCircleCategories(),
            },
            y: {label: {text: 'Average Time (ms)', position: 'outer-middle'}}
        },
        tooltip : {
            format : {
                title: function(x) { return (x + MIN_CIRCLES) + " Sets";} ,
                value: function(x) { return "" + x.toFixed(2) + " ms"; },
            },
        },
        data:precalcSpeed || getSpeedData() });

    var lastCount = 0, lastTime = 0, totalCount = 0;
    var circleCount = MIN_CIRCLES;
    paused = true;
    var element = d3.select(selector);
    var pauseIcon = element.select(".pause-button .glyphicon");


    function pause() {
        paused = true;
        element.select(".pause-label").text(" Resume");
        pauseIcon.attr("class", "glyphicon glyphicon-play");
    }

    function resume() {
        paused = false;
        element.select(".pause-label").text(" Pause");
        pauseIcon.attr("class", "glyphicon glyphicon-pause");
        testIteration();
    }

    function togglePause() {
        if (paused) {
            resume();
        } else {
            pause();
        }
    }

    var minRadius = 0.1, maxRadius = 1.0;
    function setRadii(low, hi) {
        minRadius = low;
        maxRadius = hi;
        for (var i = 0; i < algorithms.length; ++i) {
            algorithms[i].reset();
        }
        lastCount = 0;
        lastTime = 0;
        totalCount = 0;
        circleCount = MIN_CIRCLES;
        resume();
    }
    this.setRadii = setRadii;

    element.select(".pause-button").on("click",  togglePause);

    var intersectionFunction = getIntersections;
    this.intersectionFunction = function(_) {
        if (!arguments.length) return intersectionFunction;
        intersectionFunction = _;
    };

    // single iteration
    function testIteration() {
        element.select("#warning").style("display", "none");
        circleCount += 1;
        if (circleCount > MAX_CIRCLES) {
            circleCount = MIN_CIRCLES;
        }

        // slightly awkkard here, needed to work properly with deferreds,
        // returns a function that when called with the computed values
        // updates the algorithm counts.
        function update(algorithm, intersections, start) {
            return function(computed) {
                var elapsed = performance.now() - start;
                algorithm.update(intersections, computed, elapsed);
            };
        }

        var circles = generateRandomCircles(circleCount, minRadius, maxRadius);

        var intersections = intersectionFunction(circles);

        var deferred = [];
        for (var i = 0; i < algorithms.length; ++i ) {
            var updateResults = update(algorithms[i], intersections, performance.now()),
                computed = algorithms[i].method(intersections);

            if ("then" in computed) {
                deferred.push(computed.then(updateResults));
            } else {
                updateResults(computed);
            }
        }

        // wait for deferred objects to finish, possibly reload graph if
        // needed
        $.when.apply($, deferred).then(function() {
            var count = 0;
            for (var i = MIN_CIRCLES; i < MAX_CIRCLES; ++i) {
                count += algorithms[0].counts[i];
            }
            element.select("#iterations").text("Total trails: " + count);

            // constantly reloading the charts causes issues with hover
            // slow it down to every 10 seconds, with faster reloads at beginning
            if (((Math.log(count) - Math.log(lastCount))  >= 0.1) ||
                (performance.now() >= lastTime + 10 * 1000))  {
                performanceChart.load(getPerformanceData());
                speedChart.load(getSpeedData());
                lastCount = count;
                lastTime = performance.now();
            }

            if (!paused) {
                window.setTimeout(testIteration, 0);
            }
        }, function(error) {
            console.log("Failed to run performance test");
            console.log(error);
            element.select("#warning").html("<strong>Failed!</strong> Well that didn't work, is everything running ok?").style("display", "block");
            pause();
        });
    }
}
