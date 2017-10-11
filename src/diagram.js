import { select, selectAll } from "d3-selection";
import { transition } from "d3-transition";

import { venn, normalizeSolution, scaleSolution } from "./layout";
import { intersectionArea, distance, getCenter, getIntersectionPoints, containedInCircles } from "./circleintersection";
import { nelderMead } from "../node_modules/fmin/index.js";

/*global console:true*/

export function VennDiagram(option) {
	var width = 600,
		height = 350,
		padding = 15,
		duration = 1000,
		orientation = Math.PI / 2,
		normalize = true,
		wrap = true,
		styled = true,
		fontSize = null,
		orientationOrder = null,
		splitIntersection = option.splitIntersection || false,

		// mimic the behaviour of d3.scale.category10 from the previous
		// version of d3
		colourMap = {},

		// so this is the same as d3.schemeCategory10, which is only defined in d3 4.0
		// since we can support older versions of d3 as long as we don't force this,
		// I'm hackily redefining below. TODO: remove this and change to d3.schemeCategory10
		colourScheme = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"],
		colourIndex = 0,
		colours = function (key) {
			if (key in colourMap) {
				return colourMap[key];
			}
			var ret = colourMap[key] = colourScheme[colourIndex];
			colourIndex += 1;
			if (colourIndex >= colourScheme.length) {
				colourIndex = 0;
			}
			return ret;
		},
		layoutFunction = venn;

	function chart(selection) {
		var data = selection.datum();
		var solution = layoutFunction(data);
		if (normalize) {
			solution = normalizeSolution(solution,
				orientation,
				orientationOrder);
		}
		var circles = scaleSolution(solution, width, height, padding);
		var textCentres = computeTextCentres(circles, data);

		// create svg if not already existing
		selection.selectAll("svg").data([circles]).enter().append("svg");

		var svg = selection.select("svg")
			.attr("width", width)
			.attr("height", height);

		// to properly transition intersection areas, we need the
		// previous circles locations. load from elements
		var previous = {}, hasPrevious = false;
		svg.selectAll(".venn-area path").each(function (d) {
			var path = select(this).attr("d");
			if ((d.sets.length == 1) && path) {
				hasPrevious = true;
				previous[d.sets[0]] = circleFromPath(path);
			}
		});

		// interpolate intersection area paths between previous and
		// current paths
		var pathTween = function (d) {
			return function (t) {
				var c = d.sets.map(function (set) {
					var start = previous[set], end = circles[set];
					if (!start) {
						start = { x: width / 2, y: height / 2, radius: 1 };
					}
					if (!end) {
						end = { x: width / 2, y: height / 2, radius: 1 };
					}
					return {
						'x': start.x * (1 - t) + end.x * t,
						'y': start.y * (1 - t) + end.y * t,
						'radius': start.radius * (1 - t) + end.radius * t
					};
				});
				return intersectionAreaPath(c);
			};
		};

		// update data, joining on the set ids
		var nodes = svg.selectAll(".venn-area")
			.data(data, function (d) { return d.sets; });

		// create new nodes
		var enter = nodes.enter()
			.append('g')
			.attr("class", function (d) {
				return "venn-area venn-" +
					(d.sets.length == 1 ? "circle" : "intersection");
			})
			.attr("data-venn-sets", function (d) {
				return d.sets.join("_");
			});

		var enterPath = enter.append("path"),
			enterText = enter.append("text")
				.attr("class", "label")
				.text(function (d) { return label(d); })
				.attr("text-anchor", "middle")
				.attr("dy", ".35em")
				.attr("x", width / 2)
				.attr("y", height / 2);

		// get the curves of intersection
		if (splitIntersection) {
			var arcs = {},
				keys = Object.keys(circles),
				curves = {},
				max = 0;
			for (var datum = 0; datum < selection.datum().length; ++datum) {
				var each = selection.datum()[datum];
				if (each.sets.length > max) {
					max = each.sets.length;
				}
			}
			for (var lines = 2; lines <= max; lines++) {
				var tempLines = [];
				for (var dat = 0; dat < selection.datum().length; ++dat) {
					if (selection.datum()[dat].sets.length === lines) {
						tempLines.push(selection.datum()[dat]);
					}
				}
				curves[lines] = [];
				for (var tempIndex = 0; tempIndex < tempLines.length; ++tempIndex) {
					var line = tempLines[tempIndex];
					var needCircles = new Array(line.sets.length);
					for (var i = 0; i < line.sets.length; i++) {
						needCircles[i] = circles[line.sets[i]];
					}
					intersectionArea(needCircles, arcs);
					for (var arc = 0; arc < arcs.arcs.length; ++arc) {
						var temp = {};
						temp.width = arcs.arcs[arc].width;
						temp.radius = arcs.arcs[arc].circle.radius;
						temp.center = { x: arcs.arcs[arc].circle.x, y: arcs.arcs[arc].circle.y };
						temp.p1 = { x: arcs.arcs[arc].p1.x, y: arcs.arcs[arc].p1.y };
						temp.p2 = { x: arcs.arcs[arc].p2.x, y: arcs.arcs[arc].p2.y };
						temp.parentIndex = line.sets.length == 2 ? line.sets : arcs.arcs[arc].p1.parentIndex;
						curves[lines].push(temp);
					}
				}
			}
		}


		for (var len = 0; len < data.length; len++) {
			var d = data[len];
			if (d.sets.length == 1) {
				d.colour = colours(label(d));
			} else {
				d.colour = colourMap[label(d)] = colourScheme[colourIndex];
				colourIndex += 1;
				if (colourIndex >= colourScheme.length) {
					colourIndex = 0;
				}
			}
		}
		// apply minimal style if wanted
		if (styled) {
			enterPath.style("fill-opacity", "0")
				.style("fill", function (d) { return d.colour; })
				.style("fill-opacity", "0.25");

			enterText
				.style("fill", function (d) { return d.colour; });
		}

		// update existing, using pathTween if necessary
		var update = selection;
		if (hasPrevious) {
			update = selection.transition("venn").duration(duration);
			update.selectAll("path")
				.attrTween("d", pathTween);
		} else {
			update.selectAll("path")
				.attr("d", function (d) {
					if (splitIntersection && curves[2].length !== 0 && d.sets.length < max) {
						var needCurves = [],
							oppoCurves = [],
							circle = circles[d.sets],
							lines = d.sets.length === 1 ? 2 : d.sets.length;
						for (var cur = 0; cur < curves[lines].length; ++cur) {
							var each = curves[lines][cur];
							if (d.sets.length === 1) {
								if (each.parentIndex.indexOf(d.sets[0]) !== -1) {
									if (each.radius !== circle.radius) {
										needCurves.push(each);
									}
								}
							} else {
								if (each.parentIndex.indexOf(d.sets[0]) !== -1 && each.parentIndex.indexOf(d.sets[1]) !== -1) {
									needCurves.push(each);
								}
							}
						}
						if (d.sets.length > 1 && needCurves.length !== 0) {
							intersectReadjust(needCurves, curves[lines + 1]);
						}
						if (d.sets.length === 1) {
							needCurves = arcReadjust(needCurves, circle, curves, d.sets);
						}
						if (needCurves.length !== 0) {
							var order = rearrange(needCurves, circle, d.sets);
							var ret = ['\nM', needCurves[0].p1.x, needCurves[0].p1.y];
							for (var point = 0; point < order.length; ++point) {
								if (circle) {
									var r = order[point][3] === 0 ? order[point][0] : circle.radius;
									ret.push('\nA', r, r, 0, order[point][1], order[point][3], order[point][2].x, order[point][2].y);
								} else {
									ret.push('\nA', order[point][0], order[point][0], 0, order[point][1], order[point][3], order[point][2].x, order[point][2].y);
								}
							}
							return ret.join(" ");
						}
					} else
						return intersectionAreaPath(d.sets.map(function (set) { return circles[set]; }));
				});
		}

		var updateText = update.selectAll("text")
			.filter(function (d) { return d.sets in textCentres; })
			.text(function (d) { return label(d); })
			.attr("x", function (d) { return Math.floor(textCentres[d.sets].x); })
			.attr("y", function (d) { return Math.floor(textCentres[d.sets].y); });

		if (wrap) {
			if (hasPrevious) {
				// d3 4.0 uses 'on' for events on transitions,
				// but d3 3.0 used 'each' instead. switch appropiately
				if ('on' in updateText) {
					updateText.on("end", wrapText(circles, label));
				} else {
					updateText.each("end", wrapText(circles, label));
				}
			} else {
				updateText.each(wrapText(circles, label));
			}
		}

		// remove old
		var exit = nodes.exit().transition('venn').duration(duration).remove();
		exit.selectAll("path")
			.attrTween("d", pathTween);

		var exitText = exit.selectAll("text")
			.attr("x", width / 2)
			.attr("y", height / 2);

		// if we've been passed a fontSize explicitly, use it to
		// transition
		if (fontSize !== null) {
			enterText.style("font-size", "0px");
			updateText.style("font-size", fontSize);
			exitText.style("font-size", "0px");
		}


		return {
			'circles': circles,
			'textCentres': textCentres,
			'nodes': nodes,
			'enter': enter,
			'update': update,
			'exit': exit
		};
	}

	function label(d) {
		if (d.label) {
			return d.label;
		}
		if (d.sets.length == 1) {
			return '' + d.sets[0];
		}
	}

	chart.wrap = function (_) {
		if (!arguments.length) return wrap;
		wrap = _;
		return chart;
	};

	chart.width = function (_) {
		if (!arguments.length) return width;
		width = _;
		return chart;
	};

	chart.height = function (_) {
		if (!arguments.length) return height;
		height = _;
		return chart;
	};

	chart.padding = function (_) {
		if (!arguments.length) return padding;
		padding = _;
		return chart;
	};

	chart.colours = function (_) {
		if (!arguments.length) return colours;
		colours = _;
		return chart;
	};

	chart.fontSize = function (_) {
		if (!arguments.length) return fontSize;
		fontSize = _;
		return chart;
	};

	chart.duration = function (_) {
		if (!arguments.length) return duration;
		duration = _;
		return chart;
	};

	chart.layoutFunction = function (_) {
		if (!arguments.length) return layoutFunction;
		layoutFunction = _;
		return chart;
	};

	chart.normalize = function (_) {
		if (!arguments.length) return normalize;
		normalize = _;
		return chart;
	};

	chart.styled = function (_) {
		if (!arguments.length) return styled;
		styled = _;
		return chart;
	};

	chart.orientation = function (_) {
		if (!arguments.length) return orientation;
		orientation = _;
		return chart;
	};

	chart.orientationOrder = function (_) {
		if (!arguments.length) return orientationOrder;
		orientationOrder = _;
		return chart;
	};

	return chart;
}
// sometimes text doesn't fit inside the circle, if thats the case lets wrap
// the text here such that it fits
// todo: looks like this might be merged into d3 (
// https://github.com/mbostock/d3/issues/1642),
// also worth checking out is
// http://engineering.findthebest.com/wrapping-axis-labels-in-d3-js/
// this seems to be one of those things that should be easy but isn't
export function wrapText(circles, labeller) {
	return function () {
		var text = select(this),
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
			.attr("dy", function (d, i) {
				return (initial + i * lineHeight) + "em";
			});
	};
}

function circleMargin(current, interior, exterior) {
	var margin = interior[0].radius - distance(interior[0], current), i, m;
	for (i = 1; i < interior.length; ++i) {
		m = interior[i].radius - distance(interior[i], current);
		if (m <= margin) {
			margin = m;
		}
	}

	for (i = 0; i < exterior.length; ++i) {
		m = distance(exterior[i], current) - exterior[i].radius;
		if (m <= margin) {
			margin = m;
		}
	}
	return margin;
}

// compute the center of some circles by maximizing the margin of
// the center point relative to the circles (interior) after subtracting
// nearby circles (exterior)
export function computeTextCentre(interior, exterior) {
	// get an initial estimate by sampling around the interior circles
	// and taking the point with the biggest margin
	var points = [], i;
	for (i = 0; i < interior.length; ++i) {
		var c = interior[i];
		points.push({ x: c.x, y: c.y });
		points.push({ x: c.x + c.radius / 2, y: c.y });
		points.push({ x: c.x - c.radius / 2, y: c.y });
		points.push({ x: c.x, y: c.y + c.radius / 2 });
		points.push({ x: c.x, y: c.y - c.radius / 2 });
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
	var solution = nelderMead(
		function (p) { return -1 * circleMargin({ x: p[0], y: p[1] }, interior, exterior); },
		[initial.x, initial.y],
		{ maxIterations: 500, minErrorDelta: 1e-10 }).x;
	var ret = { x: solution[0], y: solution[1] };

	// check solution, fallback as needed (happens if fully overlapped
	// etc)
	var valid = true;
	for (i = 0; i < interior.length; ++i) {
		if (distance(ret, interior[i]) > interior[i].radius) {
			valid = false;
			break;
		}
	}

	for (i = 0; i < exterior.length; ++i) {
		if (distance(ret, exterior[i]) < exterior[i].radius) {
			valid = false;
			break;
		}
	}

	if (!valid) {
		if (interior.length == 1) {
			ret = { x: interior[0].x, y: interior[0].y };
		} else {
			var areaStats = {};
			intersectionArea(interior, areaStats);

			if (areaStats.arcs.length === 0) {
				ret = { 'x': 0, 'y': -1000, disjoint: true };

			} else if (areaStats.arcs.length == 1) {
				ret = {
					'x': areaStats.arcs[0].circle.x,
					'y': areaStats.arcs[0].circle.y
				};

			} else if (exterior.length) {
				// try again without other circles
				ret = computeTextCentre(interior, []);

			} else {
				// take average of all the points in the intersection
				// polygon. this should basically never happen
				// and has some issues:
				// https://github.com/benfred/venn.js/issues/48#issuecomment-146069777
				ret = getCenter(areaStats.arcs.map(function (a) { return a.p1; }));
			}
		}
	}

	return ret;
}

// given a dictionary of {setid : circle}, returns
// a dictionary of setid to list of circles that completely overlap it
function getOverlappingCircles(circles) {
	var ret = {}, circleids = [];
	for (var circleid in circles) {
		circleids.push(circleid);
		ret[circleid] = [];
	}
	for (var i = 0; i < circleids.length; i++) {
		var a = circles[circleids[i]];
		for (var j = i + 1; j < circleids.length; ++j) {
			var b = circles[circleids[j]],
				d = distance(a, b);

			if (d + b.radius <= a.radius + 1e-10) {
				ret[circleids[j]].push(circleids[i]);

			} else if (d + a.radius <= b.radius + 1e-10) {
				ret[circleids[i]].push(circleids[j]);
			}
		}
	}
	return ret;
}

export function computeTextCentres(circles, areas) {
	var ret = {}, overlapped = getOverlappingCircles(circles);
	for (var i = 0; i < areas.length; ++i) {
		var area = areas[i].sets, areaids = {}, exclude = {};
		for (var j = 0; j < area.length; ++j) {
			areaids[area[j]] = true;
			var overlaps = overlapped[area[j]];
			// keep track of any circles that overlap this area,
			// and don't consider for purposes of computing the text
			// centre
			for (var k = 0; k < overlaps.length; ++k) {
				exclude[overlaps[k]] = true;
			}
		}

		var interior = [], exterior = [];
		for (var setid in circles) {
			if (setid in areaids) {
				interior.push(circles[setid]);
			} else if (!(setid in exclude)) {
				exterior.push(circles[setid]);
			}
		}
		var centre = computeTextCentre(interior, exterior);
		ret[area] = centre;
		if (centre.disjoint && (areas[i].size > 0)) {
			console.log("WARNING: area " + area + " not represented on screen");
		}
	}
	return ret;
}

// sorts all areas in the venn diagram, so that
// a particular area is on top (relativeTo) - and
// all other areas are so that the smallest areas are on top
export function sortAreas(div, relativeTo) {

	// figure out sets that are completly overlapped by relativeTo
	var overlaps = getOverlappingCircles(div.selectAll("svg").datum());
	var exclude = {};
	for (var i = 0; i < relativeTo.sets.length; ++i) {
		var check = relativeTo.sets[i];
		for (var setid in overlaps) {
			var overlap = overlaps[setid];
			for (var j = 0; j < overlap.length; ++j) {
				if (overlap[j] == check) {
					exclude[setid] = true;
					break;
				}
			}
		}
	}

	// checks that all sets are in exclude;
	function shouldExclude(sets) {
		for (var i = 0; i < sets.length; ++i) {
			if (!(sets[i] in exclude)) {
				return false;
			}
		}
		return true;
	}

	// need to sort div's so that Z order is correct
	div.selectAll("g").sort(function (a, b) {
		// highest order set intersections first
		if (a.sets.length != b.sets.length) {
			return a.sets.length - b.sets.length;
		}

		if (a == relativeTo) {
			return shouldExclude(b.sets) ? -1 : 1;
		}
		if (b == relativeTo) {
			return shouldExclude(a.sets) ? 1 : -1;
		}

		// finally by size
		return b.size - a.size;
	});
}

export function circlePath(x, y, r) {
	var ret = [];
	ret.push("\nM", x, y);
	ret.push("\nm", -r, 0);
	ret.push("\na", r, r, 0, 1, 0, r * 2, 0);
	ret.push("\na", r, r, 0, 1, 0, -r * 2, 0);
	return ret.join(" ");
}

// inverse of the circlePath function, returns a circle object from an svg path
export function circleFromPath(path) {
	var tokens = path.split(' ');
	return {
		'x': parseFloat(tokens[1]),
		'y': parseFloat(tokens[2]),
		'radius': -parseFloat(tokens[4])
	};
}

/** returns a svg path of the intersection area of a bunch of circles */
export function intersectionAreaPath(circles) {
	var stats = {};
	intersectionArea(circles, stats);
	var arcs = stats.arcs;

	if (arcs.length === 0) {
		return "M 0 0";

	} else if (arcs.length == 1) {
		var circle = arcs[0].circle;
		return circlePath(circle.x, circle.y, circle.radius);

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
}

// check if a point is on intersection arc(s)
function containedInArcs(point, arcs) {
	for (var arc = 0; arc < arcs.length; ++arc) {
		var each = arcs[arc],
			a = each.p2.x - each.p1.x,
			b = each.p2.y - each.p1.y;
		var pointToLine = b * (point.x - each.p1.x) - a * (point.y - each.p1.y),
			centerToLine = b * (each.center.x - each.p1.x) - a * (each.center.y - each.p1.y);
		if (onSameSide([each.p1, each.p2], each.center, 0)) {
			if (pointToLine * centerToLine >= 0) {
				return false;
			}
		} else {
			if (pointToLine * centerToLine <= 0) {
				return false;
			}
		}
	}
	return true;
}

// check whether a point contain in arc(s) of the circle with center given
function containedInArc(point, arcs, center) {
	for (var arc = 0; arc < arcs.length; ++arc) {
		var each = arcs[arc],
			a = each.p2.x - each.p1.x,
			b = each.p2.y - each.p1.y;
		var pointToLine = b * (point.x - each.p1.x) - a * (point.y - each.p1.y),
			centerToLine = b * (center.x - each.p1.x) - a * (center.y - each.p1.y);
		if (onSameSide([each.p1, each.p2], center, 1)) {
			if (pointToLine * centerToLine >= 0) {
				return false;
			}
		} else {
			if (pointToLine * centerToLine <= 0) {
				return 0;
			}
		}
	}
	return true;
}

// get the order to draw
function rearrange(needCurves, center, sets) {
	var firstCurve = needCurves[0],
		wide = firstCurve.width > firstCurve.radius ? 1 : 0,
		order = [[firstCurve.radius, wide, firstCurve.p2, 0]],
		current = [firstCurve.radius, wide, firstCurve.p2, 0],
		used = [0];
	do {
		var next = nextPoint(current[2], needCurves, used);
		if (next) {
			order.push(next.slice());
			current = next;
		} else {
			var shortest = [2 * Math.PI, firstCurve.p1, firstCurve.radius];
			for (var i = 1; i < needCurves.length; i++) {
				if (used.indexOf(i) !== -1) {
					continue;
				}
				var each = needCurves[i],
					examPoint = getNearestPoint([each.p1, each.p2], current[2], center, 1),
					arcAngle = calculateRadian({ center: center, p1: current[2], p2: examPoint });
				if (!onSameSide([current[2], examPoint], center, 1)) {
					arcAngle = 2 * Math.PI - arcAngle;
				}
				if (arcAngle < shortest[0]) {
					shortest = [arcAngle, examPoint, each.radius, i];
				}
			}
			if (shortest[0] === 2 * Math.PI) {
				wide = onSameSide([current[2], firstCurve.p1], center, 1) ? 0 : 1;
				order.push([0, wide, firstCurve.p1, 1]);
				current = [0, wide, firstCurve.p1, 1];
			} else {
				wide = shortest[0] > Math.PI ? 1 : 0;
				current = [shortest[2], wide, shortest[1], 1];
				order.push(current.slice());
			}
		}
	} while (!isSamePoint(current[2], firstCurve.p1));
	return order;
}

// check the arc is big or small arc
function onSameSide(linePoints, center, sweep) {
	// get the line equation
	var first = linePoints[0],
		second = linePoints[1],
		a = second.x - first.x,
		b = second.y - first.y,
		n;
	// get the normal vector of clockwise and counter clockwise
	if (sweep === 1) {
		n = { x: center.y - first.y, y: first.x - center.x };
	} else {
		n = { x: first.y - center.y, y: center.x - first.x };
	}
	var cos = calculateRadian(undefined, [n, { x: a, y: b }]);
	if (cos <= Math.PI / 2) {
		return true;
	} else return false;
}


function calculateRadian(arc, vectos) {
	var firstVecto, secondVecto;
	if (!vectos || vectos.length === 0) {
		firstVecto = { x: arc.center.x - arc.p1.x, y: arc.center.y - arc.p1.y };
		secondVecto = { x: arc.center.x - arc.p2.x, y: arc.center.y - arc.p2.y };
	} else {
		firstVecto = vectos[0];
		secondVecto = vectos[1];
	}
	var multiple = firstVecto.x * secondVecto.x + firstVecto.y * secondVecto.y,
		cos;
	if (!vectos || vectos.length === 0) {
		cos = multiple / (distance(arc.center, arc.p1) * distance(arc.center, arc.p2));
	} else {
		cos = multiple / (Math.sqrt(vectos[0].x * vectos[0].x + vectos[0].y * vectos[0].y) * Math.sqrt(vectos[1].x * vectos[1].x + vectos[1].y * vectos[1].y));
	}
	return Math.acos(cos);
}


function nextPoint(current, arcs, used) {
	var next;
	for (var i = 0; i < arcs.length; i++) {
		if (used.indexOf(i) !== -1) {
			continue;
		}
		var each = arcs[i],
			wide = each.width > each.radius ? 1 : 0,
			sweep;
		if (isSamePoint(current, each.p1)) {
			sweep = each.sweep ? 1 : 0;
			next = [each.radius, wide, each.p2, sweep];
			used.push(i);
			return next;
		}
		if (isSamePoint(current, each.p2)) {
			sweep = each.sweep ? 1 : 0;
			next = [each.radius, wide, each.p1, sweep];
			used.push(i);
			return next;
		}
	}
}

function getNearestPoint(endPoints, start, center, sweep) {
	if (sweep === 1) {
		if (containedInArc(endPoints[1], [{ center: center, p1: start, p2: endPoints[0] }], center)) {
			return endPoints[1];
		} else return endPoints[0];
	} else {
		if (!containedInArc(endPoints[1], [{ center: center, p1: start, p2: endPoints[0] }], center)) {
			return endPoints[1];
		} else return endPoints[0];
	}
}

// readjust arcs that has intersect point
function arcReadjust(needCurves, center, curves, sets) {
	for (var i = 0; i < needCurves.length; i++) {
		var a = needCurves[i],
			circleA = { radius: a.radius, x: a.center.x, y: a.center.y };
		for (var j = 0; j < needCurves.length; j++) {
			var b = needCurves[j];
			if (i == j) {
				continue;
			}
			if (isSamePoint(a.p1, b.p2) || isSamePoint(a.p2, b.p1)) {
				continue;
			}
			var circleB = { radius: b.radius, x: b.center.x, y: b.center.y },
				intersectionPoints = getIntersectionPoints([circleA, circleB]);
			if (intersectionPoints.length == 0) {
				continue;
			}
			var needPoint = [];
			for (var point = 0; point < intersectionPoints.length; ++point) {
				if (containedInArcs(intersectionPoints[point], [a, b])) {
					needPoint.push(intersectionPoints[point]);
				}
			}
			if (needPoint.length === 0) {
				continue;
			}
			if (needPoint.length === 1) {
				var copy;
				if (containedInArc(a.p1, [b], center)) {
					copy = Object.assign({}, a);
					a.p1 = { x: needPoint[0].x, y: needPoint[0].y };
					reSize(a);
				} else if (containedInArc(a.p2, [b], center)) {
					copy = Object.assign({}, a);
					a.p2 = { x: needPoint[0].x, y: needPoint[0].y };
					reSize(a);
				}
				if (containedInArc(b.p1, [copy], center)) {
					b.p1 = { x: needPoint[0].x, y: needPoint[0].y };
					reSize(b);
				} else if (containedInArc(b.p2, [copy], center)) {
					b.p2 = { x: needPoint[0].x, y: needPoint[0].y };
					reSize(b);
				}
			}
			if (needPoint.length === 2) {
				var p1p1arc = calculateRadian({ center: center, p1: a.p1, p2: b.p1 }),
					p1p2arc = calculateRadian({ center: center, p1: a.p1, p2: b.p2 }),
					p2p1arc = calculateRadian({ center: center, p1: a.p2, p2: a.p1 }),
					p2p2arc = calculateRadian({ center: center, p1: a.p2, p2: a.p1 });
				if (!onSameSide([a.p1, b.p1], center, 1)) {
					p1p1arc = Math.PI * 2 - p1p1arc;
				}
				if (!onSameSide([a.p1, b.p2], center, 1)) {
					p1p2arc = Math.PI * 2 - p1p2arc;
				}
				if (!onSameSide([a.p2, b.p1], center, 1)) {
					p2p1arc = Math.PI * 2 - p2p1arc;
				}
				if (!onSameSide([a.p2, b.p2], center, 1)) {
					p2p2arc = Math.PI * 2 - p2p2arc;
				}
				var max = Math.max(p1p1arc, p1p2arc, p2p1arc, p2p2arc),
					replacePoint,
					newStartPoint;
				if (p1p1arc === max) {
					replacePoint = getNearestPoint(needPoint, a.p1, center, 0);
					if (isSamePoint(replacePoint, needPoint[0])) {
						newStartPoint = needPoint[1];
					} else {
						newStartPoint = needPoint[0];
					}
					if (a.width > a.radius) {
						if (onSameSide([newStartPoint, a.p2], a.center, 0)) {
							a.width = a.width - a.radius;
						}
					}
					if (b.width > b.radius) {
						if (onSameSide([newStartPoint, b.p2], b.center, 0)) {
							b.width = b.width - b.radius;
						}
					}
					curves[2].push({ width: a.width, radius: a.radius, center: a.center, p1: newStartPoint, p2: a.p2, parentIndex: a.parentIndex });
					curves[2].push({ width: b.width, radius: b.radius, center: b.center, p1: newStartPoint, p2: b.p2, parentIndex: b.parentIndex });
					b.p2 = a.p2 = { x: replacePoint.x, y: replacePoint.y };
				}
				if (p1p2arc === max) {
					replacePoint = getNearestPoint(needPoint, a.p1, center, 0);
					if (isSamePoint(replacePoint, needPoint[0])) {
						newStartPoint = needPoint[1];
					} else {
						newStartPoint = needPoint[0];
					}
					if (a.width > a.radius) {
						if (onSameSide([newStartPoint, a.p2], a.center, 0)) {
							a.width = a.width - a.radius;
						}
					}
					if (b.width > b.radius) {
						if (onSameSide([b.p1, newStartPoint], b.center, 0)) {
							b.width = b.width - b.radius;
						}
					}
					curves[2].push({ width: a.width, radius: a.radius, center: a.center, p1: newStartPoint, p2: a.p2, parentIndex: a.parentIndex });
					curves[2].push({ width: b.width, radius: b.radius, center: b.center, p1: b.p1, p2: newStartPoint, parentIndex: b.parentIndex });
					b.p1 = a.p2 = { x: replacePoint.x, y: replacePoint.y };
				}
				if (p2p1arc === max) {
					replacePoint = getNearestPoint(needPoint, a.p2, center, 0);
					if (isSamePoint(replacePoint, needPoint[0])) {
						newStartPoint = needPoint[1];
					} else {
						newStartPoint = needPoint[0];
					}
					if (a.width > a.radius) {
						if (onSameSide([a.p1, newStartPoint], a.center, 0)) {
							a.width = a.width - a.radius;
						}
					}
					if (b.width > b.radius) {
						if (onSameSide([newStartPoint, b.p2], b.center, 0)) {
							b.width = b.width - b.radius;
						}
					}
					curves[2].push({ width: a.width, radius: a.radius, center: a.center, p1: a.p1, p2: newStartPoint, parentIndex: a.parentIndex });
					curves[2].push({ width: b.width, radius: b.radius, center: b.center, p1: newStartPoint, p2: b.p2, parentIndex: b.parentIndex });
					b.p2 = a.p1 = { x: replacePoint.x, y: replacePoint.y };
				}
				if (p2p2arc === max) {
					replacePoint = getNearestPoint(needPoint, a.p2, center, 0);
					if (isSamePoint(replacePoint, needPoint[0])) {
						newStartPoint = needPoint[1];
					} else {
						newStartPoint = needPoint[0];
					}
					if (a.width > a.radius) {
						if (onSameSide([a.p1, newStartPoint], a.center, 0)) {
							a.width = a.width - a.radius;
						}
					}
					if (b.width > b.radius) {
						if (onSameSide([b.p1, newStartPoint], b.center, 0)) {
							b.width = b.width - b.radius;
						}
					}
					curves[2].push({ width: a.width, radius: a.radius, center: a.center, p1: a.p1, p2: newStartPoint, parentIndex: a.parentIndex });
					curves[2].push({ width: b.width, radius: b.radius, center: b.center, p1: b.p1, p2: newStartPoint, parentIndex: b.parentIndex });
					b.p1 = a.p1 = { x: replacePoint.x, y: replacePoint.y };
				}
			}
		}
	}
	return needCurves;
}

function isSamePoint(p1, p2) {
	if (Math.abs(p1.x - p2.x) < 1e-10 && Math.abs(p1.y - p2.y) < 1e-10) {
		return true;
	} else return false;
}

function reSize(arc) {
	if (arc.width > arc.radius) {
		if (onSameSide([arc.p1, arc.p2], arc.center, 0)) {
			arc.width = arc.width - arc.radius;
		}
	}
}

// get the remain arc missing of the intersrction
function intersectReadjust(needCurves, curves) {
	// get the missing part starting and ending point
	// by going through all the curves
	// till there no more curve connect to 
	// the current start and end point
	var examPoint = [],
		start = needCurves[0].p1,
		end = needCurves[0].p2,
		used = [0],
		nextstart;
	while (used.length !== needCurves.length) {
		nextstart = nextPoint(start, needCurves, used);
		if (nextstart) {
			start = nextstart[2];
		}
		var next = nextPoint(end, needCurves, used);
		if (next) {
			end = next[2];
		}
		if (!nextstart && !next) {
			break;
		}
	}

	// in case the intersection cut by 2 other cuves
	// there will be curves that not connect to those found before
	if (used.length !== needCurves.length) {
		start = [start];
		end = [end];
		var unUsed = [];
		for (var i = 0; i < needCurves.length; i++) {
			if (used.indexOf(i) === -1) {
				unUsed.push(needCurves[i]);
			}
		}
		used = [0];
		var otherStart = unUsed[0].p1,
			otherEnd = unUsed[0].p2,
			nextend;
		while (used.length !== unUsed) {
			nextstart = nextPoint(otherStart, unUsed, used);
			if (nextstart) {
				otherStart = nextstart[2];
			}
			nextend = nextPoint(otherEnd, unUsed, used);
			if (nextend) {
				otherEnd = nextend[2];
			}
			if (!nextend && !nextstart) {
				break;
			}
		}
		start.push(otherStart);
		end.push(otherEnd);
	}

	// search for the missing arc(s) of the intersection if the intersection has not compvare
	if (!isSamePoint(needCurves[0].p1, needCurves[1].p2) || !isSamePoint(needCurves[0].p2, needCurves[1].p1)) {
		var each, curv;
		if (start.constructor !== Array) {
			for (curv = 0; curv < curves.length; ++curv) {
				each = curves[curv];
				if (isSamePoint(each.p1, start)) {
					if (isSamePoint(each.p2, end)) {
						each.sweep = true;
						needCurves.push(each);
						break;
					}
				}
				if (isSamePoint(each.p2, start)) {
					if (isSamePoint(each.p1, end)) {
						each.sweep = true;
						needCurves.push(each);
						break;
					}
				}
			}
		} else {
			for (curv = 0; curv < curves.length; ++curv) {
				each = curves[curv];
				if (isSamePoint(each.p1, start[0])) {
					if (isSamePoint(each.p2, end[1])) {
						each.sweep = true;
						needCurves.push(each);
						continue;
					}
				}
				if (isSamePoint(each.p2, start[0])) {
					if (isSamePoint(each.p1, end[1])) {
						each.sweep = true;
						needCurves.push(each);
						continue;
					}
				}
				if (isSamePoint(each.p1, start[1])) {
					if (isSamePoint(each.p2, end[0])) {
						each.sweep = true;
						needCurves.push(each);
						continue;
					}
				}
				if (isSamePoint(each.p2, start[1])) {
					if (isSamePoint(each.p1, end[0])) {
						each.sweep = true;
						needCurves.push(each);
						continue;
					}
				}
			}
		}
	}
}