var tape = require("tape"),
    venn = require("../");

var SMALL = 1e-5;

function nearlyEqual(test, left, right, tolerance, message) {
    message = message || "nearlyEqual";
    tolerance = tolerance || SMALL;
    test.ok(Math.abs(left - right) < tolerance,
            message + ": " + left + " ~== " + right);
}

function lessThan(test, left, right, message) {
    message = message || "lessThan";
    test.ok(left < right, message + ": " + left + " < " + right);
}

tape("greedyLayout", function(test) {
    var areas = [{sets: [0], size:0.7746543297103429},
                 {sets: [1], size:0.1311252856844238},
                 {sets: [2], size:0.2659942131443344},
                 {sets: [3], size:0.44600866168641723},
                 {sets: [0,1], size:0.02051532092950205},
                 {sets: [0,2], size:0},
                 {sets: [0,3], size:0},
                 {sets: [1,2], size:0},
                 {sets: [1,3], size:0.07597023820511245},
                 {sets: [2,3], size:0}];

    var circles = venn.greedyLayout(areas)
        loss = venn.lossFunction(circles, areas);
    nearlyEqual(test, loss, 0);


    areas = [{sets: [0], size: 0.5299368855059736},
             {sets: [1], size: 0.03364187025606481},
             {sets: [2], size: 0.3121450394871512},
             {sets: [3], size: 0.0514397361783036},
             {sets: [0,1], size: 0.013912447645582351},
             {sets: [0,2], size: 0.005903647141469598},
             {sets: [0,3], size: 0.0514397361783036},
             {sets: [1,2], size: 0.012138157839477597},
             {sets: [1,3], size: 0.008010688232481479},
             {sets: [2,3], size: 0}];

    circles = venn.greedyLayout(areas);
    loss = venn.lossFunction(circles, areas);
    nearlyEqual(test, loss, 0);

    // one small circle completely overlapped in the intersection
    // area of two larger circles
    areas= [{sets: [0], size: 1.7288584050841396},
           {sets: [1], size: 0.040875831658950056},
           {sets: [2], size: 2.587146019782323},
           {sets: [0,1], size: 0.040875831658950056},
           {sets: [0,2], size: 0.5114617575187569},
           {sets: [1,2], size: 0.040875831658950056}];

    circles = venn.greedyLayout(areas);
    loss = venn.lossFunction(circles, areas);
    nearlyEqual(test, loss, 0);
    test.end();
});


tape("fmin", function(test) {
    // minimize simple 1 diminesial quadratic
    var loss = function(values) { return (values[0] - 10) * (values[0] - 10); }
    var solution = venn.fmin(loss, [0], {minErrorDelta:1e-10}).solution;
    nearlyEqual(test, solution[0], 10, 1e-10);
    test.end();
});

tape("minimizeConjugateGradient", function(test) {
    // minimize simple 1 diminesial quadratic
    var loss = function(x, xprime) { 
        xprime[0] = 2 * (x[0] - 10);
        return (x[0] - 10) * (x[0] - 10); 
    }
    var solution = venn.minimizeConjugateGradient(loss, [0]).x;
    nearlyEqual(test, solution[0], 10, 1e-10);
    test.end();
});

tape("circleIntegral", function(test) {
    nearlyEqual(test, venn.circleIntegral(10, 0), 0, SMALL,
        "empty circle test");
    nearlyEqual(test, venn.circleIntegral(10, 10),  Math.PI * 10 * 10 / 2,
        SMALL, "half circle test");
    test.end();
});

tape("circleArea", function(test) {
    nearlyEqual(test, venn.circleArea(10,0), 0, SMALL, "empty circle test");
    nearlyEqual(test, venn.circleArea(10, 10), Math.PI*10*10/2, SMALL,
        "half circle test");
    nearlyEqual(test, venn.circleArea(10, 20), Math.PI*10*10, SMALL,
        "full circle test");
    test.end();
});

tape("circleOverlap", function(test) {
    nearlyEqual(test, venn.circleOverlap(10, 10, 200), 0, SMALL,
        "nonoverlapping circles test");

    nearlyEqual(test, venn.circleOverlap(10, 10, 0), Math.PI*10*10, SMALL,
        "full overlapping circles test");

    nearlyEqual(test, venn.circleOverlap(10, 5, 5), Math.PI * 5 * 5, SMALL);
    test.end();
});

tape("distanceFromIntersectArea", function(test) {
    function testDistanceFromIntersectArea(r1, r2, overlap) {
        var distance = venn.distanceFromIntersectArea(r1, r2,
                                                      overlap);
        nearlyEqual(test, venn.circleOverlap(r1, r2, distance),
                     overlap,
                     1e-2);
    }

    testDistanceFromIntersectArea(1.9544100476116797,
                                  2.256758334191025,
                                  11);

    testDistanceFromIntersectArea(111.06512962798197,
                                  113.32348546565727,
                                  1218);

    testDistanceFromIntersectArea(44.456564007075,
                                  149.4335753619362,
                                  2799);

    testDistanceFromIntersectArea(592.890,
                                  134.75,
                                  56995);

    testDistanceFromIntersectArea(139.50778247443944,
                                  32.892784970851956,
                                  3399);

    testDistanceFromIntersectArea(4.886025119029199,
                                  5.077706251929807,
                                  75);
    test.end();
});

tape("circleCircleIntersection", function(test) {
    var testIntersection = function(p1, p2, msg) {
        var points = venn.circleCircleIntersection(p1, p2);
        // make sure that points are appropiately spaced
        for (var i = 0; i < points.length; i++) {
            var point = points[i];
            nearlyEqual(test, venn.distance(point, p1),
                        p1.radius, SMALL,
                        msg + ": test distance to p1 for point " + i);
            nearlyEqual(test, venn.distance(point, p2),
                        p2.radius, SMALL,
                        msg + ": test distance to p2 for point " + i );
        }

        return points;
    }

    // fully contained
    test.equal(venn.circleCircleIntersection({x:0, y:3, radius:10},
                                             {x:3, y:0, radius:20}).length,
               0, "fully contained test");

    // fully disjoint
    test.equal(venn.circleCircleIntersection({x:0, y:0, radius:10},
                                             {x:21, y:0, radius:10}).length,
               0, "fully disjoint test");

    // midway between 2 points on y axis
    var points = testIntersection({x:0, y:0, radius:10},
                                  {x:10, y:0, radius:10},
                                  "test midway intersection");
    test.equal(points.length, 2);
    nearlyEqual(test, points[0].x, 5);
    nearlyEqual(test, points[1].x, 5);
    nearlyEqual(test, points[0].y, -1 * points[1].y);

    // failing case from input
    var points = testIntersection({radius: 10, x: 15, y: 5},
                                  {radius: 10, x: 20, y: 0},
                                  "test intersection2");
    test.equal(points.length, 2);
    test.end();
});

tape("disjointCircles", function(test) {
    // each one of these circles overlaps all the others, but the total overlap is still 0
    var circles =  [{"x":0.909, "y":0.905, "radius":0.548},
                    {"x":0.765, "y":0.382, "radius":0.703},
                    {"x":0.630, "y":0.019, "radius":0.449},
                    {"x":0.210, "y":0.755, "radius":0.656},
                    {"x":0.276, "y":0.723, "radius":1.145},
                    {"x":0.141, "y":0.585, "radius":0.419}];

    var area = venn.intersectionArea(circles);
    test.ok(area == 0);

    // no intersection points, but the smallest circle is completely overlapped by each of the others
    circles = [{"x":0.426, "y":0.882,"radius":0.944},
               {"x":0.240, "y":0.685,"radius":0.992},
               {"x":0.010, "y":0.909,"radius":1.161},
               {"x":0.540, "y":0.475,"radius":0.410}];

    test.ok(circles[3].radius * circles[3].radius * Math.PI == venn.intersectionArea(circles));
    test.end();
});


tape("randomFailures", function(test) {
    var circles = [{"x":0.501,"y":0.320,"radius":0.629},
                   {"x":0.945,"y":0.022,"radius":1.015},
                   {"x":0.021,"y":0.863,"radius":0.261},
                   {"x":0.528,"y":0.090,"radius":0.676}],
        area = venn.intersectionArea(circles);

    test.ok(Math.abs(area - 0.0008914) < 0.0001, area);
    test.end();
});

tape("computeTextCentre", function(test) {
    var center = venn.computeTextCentre([{x: 0, y: 0, radius: 1}], []);
    nearlyEqual(test, center.x, 0);
    nearlyEqual(test, center.y, 0);

    var center = venn.computeTextCentre([{x: 0, y: 0, radius: 1}],  [{x: 0, y: 1, radius: 1}]);
    nearlyEqual(test, center.x, 0, 1e-4);
    nearlyEqual(test, center.y, -.5);
    test.end();
});

tape("normalizeSolution", function(test) {
    // test two circles that are far apart
    var solution = [{x: 0, y: 0, radius: .5}, 
                    {x: 1e10, y:0, radius: 1.5}];
    
    // should be placed close together
    var normalized = venn.normalizeSolution(solution);
    // distance should be 2, but we space things out
    lessThan(test, venn.distance(normalized[0], normalized[1]), 2.1);
    test.end();
});

tape("disjointClusters", function(test) {
    var input = [
        {
            "x": 0.8047033110633492,
            "y": 0.9396705999970436,
            "radius": 0.47156485118903224,
        },
        {
            "x": 0.7961132447235286,
            "y": 0.014027722179889679,
            "radius": 0.14554832570720466,
        },
        {
            "x": 0.28841276094317436,
            "y": 0.98081015329808,
            "radius": 0.9851036085514352,
        },
        {
            "x": 0.7689983483869582,
            "y": 0.2899463507346809,
            "radius": 0.7210563338827342,
        }
    ];

    var clusters = venn.disjointCluster(input);
    test.equal(clusters.length, 1);
    test.end();
});
