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
You can [view this example here](http://benfred.github.io/venn.js/examples/simple.html)

##### Dynamic layout

To have a layout that reacts to a change in input, you just need to recompute the areas and call updateD3Diagram to do the transition:

```javascript
// draw the initial set
var sets = venn.venn(getSets(), getSetIntersections());
venn.drawD3Diagram(d3.select(".dynamic"), sets, w, h);

// redraw the sets on any change in input
d3.selectAll("input").on("change", function() {
    var sets = venn.venn(getSets(), getSetIntersections());
    venn.updateD3Diagram(d3.select(".dynamic"), sets);
});
```

You can [view this example here](http://benfred.github.io/venn.js/examples/dynamic.html)

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
You can [view this example here](http://benfred.github.io/venn.js/examples/mds.html)

Released under the MIT License.
