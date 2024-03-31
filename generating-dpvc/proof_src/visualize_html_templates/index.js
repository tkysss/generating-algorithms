// Copyright 2021 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/bar-chart
function BarChart(data, {
  x = (d, i) => i, // given d in data, returns the (ordinal) x-value
  y = d => d, // given d in data, returns the (quantitative) y-value
  title, // given d in data, returns the title text
  marginTop = 30, // the top margin, in pixels
  marginRight = 0, // the right margin, in pixels
  marginBottom = 20, // the bottom margin, in pixels
  marginLeft = 60, // the left margin, in pixels
  width, // the outer width of the chart, in pixels
  height = 400, // the outer height of the chart, in pixels
  xDomain, // an array of (ordinal) x-values
  xRange, // [left, right]
  yType = d3.scaleLinear, // y-scale type
  yDomain, // [ymin, ymax]
  yRange = [height - marginBottom, marginTop], // [bottom, top]
  xPadding = 0.1, // amount of x-range to reserve to separate bars
  yFormat, // a format specifier string for the y-axis
  yLabel, // a label for the y-axis
  color = "currentColor" // bar fill color
} = {}) {
  // Compute values.
  const X = d3.map(data, x);
  const Y = d3.map(data, y);

  if(width === undefined) {
    width = X.length * 30 + marginLeft + marginRight;
  }
  if(xRange === undefined) {
    xRange = [marginLeft, width - marginRight]
  }

  // Compute default domains, and unique the x-domain.
  if (xDomain === undefined) xDomain = X;
  if (yDomain === undefined) yDomain = [0, d3.max(Y)];
  xDomain = new d3.InternSet(xDomain);

  // Omit any data not present in the x-domain.
  const I = d3.range(X.length).filter(i => xDomain.has(X[i]));

  // Construct scales, axes, and formats.
  const xScale = d3.scaleBand(xDomain, xRange).padding(xPadding);
  const yScale = yType(yDomain, yRange);
  const xAxis = d3.axisBottom(xScale).tickSizeOuter(0);
  const yAxis = d3.axisLeft(yScale).ticks(height / 40, yFormat);

  // Compute titles.
  if (title === undefined) {
    const formatValue = yScale.tickFormat(100, yFormat);
    title = i => `${X[i]}\n${formatValue(Y[i])}`;
  } else {
    const O = d3.map(data, d => d);
    const T = title;
    title = i => T(O[i], i, data);
  }

  const svg = d3.create("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("viewBox", [0, 0, width, height])
      .attr("style", "height: auto; height: intrinsic;");

  svg.append("g")
      .attr("transform", `translate(${marginLeft},0)`)
      .call(yAxis)
      .call(g => g.select(".domain").remove())
      .call(g => g.selectAll(".tick line").clone()
          .attr("x2", width - marginLeft - marginRight)
          .attr("stroke-opacity", 0.1))
      .call(g => g.append("text")
          .attr("x", -marginLeft)
          .attr("y", 10)
          .attr("fill", "currentColor")
          .attr("text-anchor", "start")
          .text(yLabel));

  const bar = svg.append("g")
      .attr("fill", color)
    .selectAll("rect")
    .data(I)
    .join("rect")
      .attr("x", i => xScale(X[i]))
      .attr("y", i => yScale(Y[i]))
      .attr("height", i => yScale(0) - yScale(Y[i]))
      .attr("width", xScale.bandwidth());

  if (title) {
    bar.append("title")
      .text(title);

const barLabels = svg.append("g")
    .selectAll('text')
    .data(I)
    .join('text')
    .text(i => Y[i])
      .attr("fill", "currentColor")
      .attr("text-anchor", "middle")
      .attr("font-size", 10)
      .attr("font-family", 'sans-serif')
      .attr("x", i => xScale(X[i]) + xScale.bandwidth() / 2)
      .attr("y", i => yScale(Y[i]) - 2)


  }

  svg.append("g")
      .attr("transform", `translate(0,${height - marginBottom})`)
      .call(xAxis);

  return svg.node();
}



function init() {
    d3.select('#rule_tree_depth_histogram').node().appendChild(BarChart(
        Object.entries(_global['algorithm_metadata']['rule_tree_depth_to_cnt']).map(([k, v]) => ({key: k, value: v})),
        {
            x: d => d.key,
            y: d => d.value,
            height: 100,
            color: "steelblue"
    }))

    d3.select('#rule_walk_length_histogram').node().appendChild(BarChart(
        Object.entries(_global['algorithm_metadata']['rule_walk_length_to_cnt']).map(([k, v]) => ({key: k, value: v})),
        {
            x: d => d.key,
            y: d => d.value,
            height: 100,
            color: "steelblue"
    }))

    d3.select('#rule_tree_depth_rule_walk_length_histogram').node().appendChild(BarChart(
        Object.entries(_global['algorithm_metadata']['rule_tree_depth_rule_walk_length_to_cnt']).map(([k, v]) => ({key: k, value: v})),
        {
            x: d => d.key,
            y: d => d.value,
            height: 100,
            color: "steelblue"
    }))

    d3.select('#rule_size_histogram').node().appendChild(BarChart(
        Object.entries(_global['algorithm_metadata']['rule_size_to_cnt']).map(([k, v]) => ({key: k, value: v})),
        {
            x: d => d.key,
            y: d => d.value,
            height: 100,
            color: "steelblue"
    }))
    d3.select('#branching_factor_histogram').node().appendChild(BarChart(
        Object.entries(_global['algorithm_metadata']['branching_rule_bf_rounded_hist_to_cnt']).map(([k, v]) => ({key: k, value: v})),
        {
            x: d => d.key,
            y: d => d.value,
            height: 150,
            color: "steelblue"
    }))




    console.log(_global)
}



init();
