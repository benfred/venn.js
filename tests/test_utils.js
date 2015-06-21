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
