venn.js
=======

A javascript library for laying out area proportional venn and euler diagrams.

Details of how this library works can be found on the [blog
post](http://www.benfrederickson.com/2013/05/09/venn-diagrams-with-d3.js.html)
I wrote about this. There are also more examples on that page.

#### Usage

This library depends on [d3.js](http://d3js.org/) to display the venn
diagrams.

##### Simple layout

To lay out a simple diagram, just define the sets and their sizes along with the sizes 
of all the set intersection. Calling 'venn.venn' will position the sets such
that the areas of each region are proportional to the sizes, and
'venn.drawD3Diagram' will display this diagram:

```javascript
// define sets and set set intersections
var sets = [{label: "A", size: 10}, {label: "B", size: 10}],
    overlaps = [{sets: [0,1], size: 2}];

// get positions for each set
sets = venn.venn(sets, overlaps);

// draw the diagram in the 'simple_example' div
venn.drawD3Diagram(d3.select(".simple_example"), sets, 300, 300);
```
[View this example ](http://benfred.github.io/venn.js/examples/simple.html)

##### Changing the Style

To change the style of the venn diagram, use D3 to set attributes on the 'text' and 'circle' objects that are returned from the drawD3Diagram call:

```javascript
var diagram = venn.drawD3Diagram(d3.select("#rings"),
                                 venn.venn(sets, overlaps), 
                                 500, 500);

// change the colours, add a thick border, remove fill
var colours = ['black', 'red', 'blue', 'green']
diagram.circles.style("fill-opacity", 0)
               .style("stroke-width", 10)
               .style("stroke-opacity", .5)
                .style("stroke", function(d,i) { return colours[i]; });

// make the font big and light
diagram.text.style("fill", function(d,i) { return colours[i]})
            .style("font-size", "24px")
            .style("font-weight", "100");
```

[View this example, along with other possible styles](http://benfred.github.io/venn.js/examples/styled.html)

##### Dynamic layout

To have a layout that reacts to a change in input, you just need to recompute the areas and call updateD3Diagram to do the transition:

```javascript
// draw the initial set
var sets = venn.venn(getSets(), getSetIntersections());
var diagram = venn.drawD3Diagram(d3.select(".dynamic"), sets, w, h);

// redraw the sets on any change in input
d3.selectAll("input").on("change", function() {
    var sets = venn.venn(getSets(), getSetIntersections());
    venn.updateD3Diagram(diagram, sets);
});
```

[View this example](http://benfred.github.io/venn.js/examples/dynamic.html)

##### Making the diagram interactive

Making the diagram interactive is basically the same idea as changing the style: just add event listeners to the returned elements. To change the text size and circle colours on mouseover:

```javascript
diagram.nodes
    .on("mouseover", function(d, i) {
        var node = d3.select(this).transition();
        node.select("circle").style("fill-opacity", .1);
        node.select("text").style("font-size", "36px");
    })
    .on("mouseout", function(d, i) {
        var node = d3.select(this).transition();
        node.select("circle").style("fill-opacity", 0);
        node.select("text").style("font-size", "24px");
    });
```
[View this example](http://benfred.github.io/venn.js/examples/interactive.html)

##### Adding tooltips to the intersection areas

The intersection areas aren't drawn by default, but there is some code
included here to draw a svg path element around the intersection areas. To add
a tooltip to the intersection area, use the 'intersectionAreaPath' method to
define the area, and then add appropiate events to handle the tooltip:

```javascript
diagram.svg.select("g").selectAll("path")
    .data(overlaps)
    .enter()
    .append("path")
    .attr("d", function(d) { 
        return venn.intersectionAreaPath(d.sets.map(function(j) { return sets[j]; })); 
    })
    .style("fill-opacity","0")
    .style("fill", "black")
    .style("stroke-opacity", 0)
    .style("stroke", "white")
    .style("stroke-width", "2")
    .on("mouseover", function(d, i) {
        d3.select(this).transition()
            .style("fill-opacity", .1)
            .style("stroke-opacity", 1);
        tooltip.transition().style("opacity", .9);
        tooltip.text(d.size + " users");
    })
    .on("mouseout", function(d, i) {
        d3.select(this).transition()
            .style("fill-opacity", 0)
            .style("stroke-opacity", 0);
        tooltip.transition().style("opacity", 0);
    })
    .on("mousemove", function() {
        tooltip.style("left", (d3.event.pageX) + "px")
               .style("top", (d3.event.pageY - 28) + "px");
    })
```
[View this example](http://benfred.github.io/venn.js/examples/intersection_tooltip.html)

##### MDS Layout

In most cases the greedy initial layout does a good job of positioning the
sets, but there are cases where it breaks down. One case is detailed in [this
blog post](http://www.benfrederickson.com/2013/05/16/multidimensional-scaling.html),
and it can be better laid out using [multidimensional
scaling](https://en.wikipedia.org/wiki/Multidimensional_scaling) to generate
the initial layout.

To enable this just include the [mds.js](http://github.com/benfred/mds.js)
and [numeric.js](http://numericjs.com) libraries first, and then generate the venn positions by calling:

```javascript
sets = venn.venn(sets, overlaps, {layoutFunction: venn.classicMDSLayout});
```
[View this example](http://benfred.github.io/venn.js/examples/mds.html)

Released under the MIT License.
