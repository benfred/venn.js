venn.js
=======

A javascript library for laying out area proportional venn and euler diagrams.


Usage
-----

This library depends on [d3.js](http://d3js.org/) to display the venn
diagrams.

##### Simple layout

To lay out a simple diagram, just define the sets and their sizes along with
with the sizes of all the set intersection. Calling 'venn.venn' will position the sets such that the areas of each region are proportional to the sizes, and venn.drawD3Diagram will display this diagram:

```javascript
// define sets and set set intersections
var sets = [{label: "A", size: 10}, {label: "B", size: 10}],
    overlaps = [{sets: [0,1], size: 2}];

// get positions for each set
sets = venn.venn(sets, overlaps);

// draw the diagram in the 'simple_example' div
venn.drawD3Diagram(d3.select(".simple_example"), sets, 300, 300);
```
You can [view this example here](http://raw.github.com/benfred/venn.js/master/examples/simple.js)

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

You can [view this example here](http://raw.github.com/benfred/venn.js/master/examples/dynamic.js)

Other examples of this libray in action can be found on the [blog
post](http://www.benfrederickson.com/2013/05/09/venn-diagrams-with-d3.js.html)
I wrote about this. That blog post also contains details on the approach taken
to compute the best layout for each venn diagram.

Released under the MIT License.
