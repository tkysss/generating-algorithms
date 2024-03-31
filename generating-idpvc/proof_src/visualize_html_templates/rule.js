const isoColors = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"
]
const redColor = '#dc3545';
const blueColor = '#0d6efd';

function drawGraphImage(graph, containerId, nodeColors) {
  // const dim = d3.select(containerId).node().getBoundingClientRect().width

  const width = 400;
  const height = 400;

  const nodes = d3.range(graph.n).map(u=>({id: u}));
  const links = graph.edges.map(([u, v])=>({source: u, target: v}));

  const forceNode = d3.forceManyBody().strength(-250);
  const forceLink = d3.forceLink(links).strength(0.8);

  const nodeRadius = 15;

  const simulation = d3.forceSimulation(nodes)
      .force("link", forceLink)
      .force("charge", forceNode)
      .force("center",  d3.forceCenter(width / 2, height / 2))
      .on("tick", ticked);

  const svg = d3.select(containerId).append("svg")
      .attr("width", width)
      .attr("height", height)
      // .attr("viewBox", [-width / 2, -height / 2, width, height]);

  const link = svg.append("g")
    .attr("stroke", 'black')
    .attr("stroke-width", 1)
    .selectAll("line")
    .data(links)
    .join("line");

  const node = svg.append("g")
    .selectAll("g")
    .data(nodes)
    .enter().append("g")
    .call(drag(simulation));

  const nodeCircle = node.append("circle")
    .attr("r", nodeRadius)
    .attr("stroke", 'black')
    .attr("stroke-width", 2)
    .attr('fill', d => nodeColors[d.id]);

  const nodeLabels = node.append("text")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central")
    .text(d => d.id);

  function ticked() {
    function clipx(x) {
      return Math.max(nodeRadius*2, Math.min(width - nodeRadius*2, x));
    }
    function clipy(y) {
      return Math.max(nodeRadius*2, Math.min(height - nodeRadius*2, y));
    }
    node
      .attr("transform", function(d) {
          return "translate(" + clipx(d.x) + "," + clipy(d.y) + ")";
        })

    link
      .attr("x1", d => clipx(d.source.x))
      .attr("y1", d => clipy(d.source.y))
      .attr("x2", d => clipx(d.target.x))
      .attr("y2", d => clipy(d.target.y));
  }

  function drag(simulation) {
    function dragstarted(event) {
      if (!event.active) simulation.alphaTarget(0.3).restart();
      event.subject.fx = event.subject.x;
      event.subject.fy = event.subject.y;
    }

    function dragged(event) {
      event.subject.fx = event.x;
      event.subject.fy = event.y;
    }

    function dragended(event) {
      if (!event.active) simulation.alphaTarget(0);
      event.subject.fx = null;
      event.subject.fy = null;
    }

    return d3.drag()
      .on("start", dragstarted)
      .on("drag", dragged)
      .on("end", dragended);
  }

}

_renderModalCache = {};

function renderSubsetDetailModal(data) {
  drawGraphWithRedVertices(_global.graph, '#subsetOrigGraphImage_' + data['subset_idx']);

  if(data.subset.type === 'not_solution') {
    const nodeColors = {}
    d3.range(_global.graph.n).forEach(v => {nodeColors[v] = '#aaa'});
    data.subset.subset.forEach(v => {nodeColors[v] = '#eee'})
    data.subset.path.forEach(v => {nodeColors[v] = 'yellow'})

    drawGraphImage(_global.graph, '#subsetPathGraphImage_' + data['subset_idx'], nodeColors);
  }
  else if(data.subset.type === 'solution_but_dominated') {
    {
    const nodeColors = {}
    d3.range(_global.graph.n).forEach(v => {nodeColors[v] = '#aaa'});
    data.subset.subset.forEach(v => {nodeColors[v] = '#eee'})
    data.subset.red_vertices_subset.forEach(v => {nodeColors[v] = redColor})
    d3.intersection(data.subset.subset, data.subset.red_vertices_subset).forEach(v => {nodeColors[v] = 'purple'})
    drawGraphImage(_global.graph, '#subsetRVSubsetIntesectionGraphImage_' + data['subset_idx'], nodeColors);
    }

    {
    const nodeColors = {}
    d3.range(_global.graph.n).forEach(v => {nodeColors[v] = '#aaa'});
    data.subset.subset.forEach(v => {nodeColors[v] = '#eee'})
    data.subset.red_vertices_subset.forEach(v => {nodeColors[v] = redColor})
    d3.difference(data.subset.dominated_by, data.subset.subset).forEach(v => {nodeColors[v] = 'green'})
    drawGraphImage(_global.graph, '#subsetDominationGraphImage_' + data['subset_idx'], nodeColors);
    }
  }
}

function renderExpansionDetailModal(data) {
  drawGraphWithRedVertices(_global.graph, '#expansionOrigGraphImage_' + data['expansion_idx']);
  {
    const nodeColors = getGraphRedVerticesNodeColors(_global.graph);
    nodeColors[_global.graph.n] = 'yellow';
    drawGraphImage(data.expansion.expansion_graph, '#expansionExpandedGraphImage_' + data['expansion_idx'], nodeColors);
  }

  if(data.expansion.expansion_result === 'eliminated' || data.expansion.expansion_result === 'next') {

    const isoColors = d3['schemeSet1'].concat(d3['schemeSet2'], d3['schemeSet3']);
    const expandedColors = {}
    const targetColors = {}

    const witnessing_isomorphism = data.expansion.result_data.witnessing_isomorphism;
    const target_graph = data.expansion.result_data.target_graph;

    Object.entries(witnessing_isomorphism).forEach(([e_v, t_v], idx)=>{
      expandedColors[+e_v] = isoColors[idx];
      targetColors[+t_v] = isoColors[idx];
    });

    {
      const nodeColors = {}
      d3.range(data.expansion.expansion_graph.n).forEach(v => {nodeColors[v] = '#eee'});
      Object.entries(expandedColors).forEach(([v, color])=>nodeColors[v] = color);
      drawGraphImage(data.expansion.expansion_graph, '#expansionExpandedIsomorphismGraphImage_' + data['expansion_idx'], nodeColors);
    }

    {
      const nodeColors = {}
      d3.range(target_graph.n).forEach(v => {nodeColors[v] = '#eee'});
      Object.entries(targetColors).forEach(([v, color])=>nodeColors[v] = color);
      drawGraphImage(target_graph, '#expansionTargetIsomorphismGraphImage_' + data['expansion_idx'], nodeColors);
    }

  }
}

function renderModal(data) {
  if(_renderModalCache[data.render_id]) return;
  _renderModalCache[data.render_id]=true;

  if(data.type == 'subset') {
    renderSubsetDetailModal(data);
  }
  else if(data.type == 'expansion') {
    renderExpansionDetailModal(data);
  }
  else {
    console.log('renderModal', data);
  }
}

function getGraphRedVerticesNodeColors(graph) {
  return Object.fromEntries(
    d3.range(graph.n).map(v => [v, ((graph.red_vertices || []).indexOf(v) != -1) ? redColor: blueColor])
  );
}

function drawGraphWithRedVertices(graph, containerId) {
  drawGraphImage(graph, containerId, getGraphRedVerticesNodeColors(graph));
}

function renderReductionRuleRedStarReductionMainGraphImage() {
    const node = document.getElementById('reductionRuleRedStarReductionMainGraphImage');
    if(!node) return;
    const nodeColors = {};
    d3.range(_global.graph.n).forEach(v => {nodeColors[v] = '#aaa'});
    _global.rule.kcenter_vertices.forEach(v => {nodeColors[v] = 'purple'});
    _global.rule.red_rays_vertices.forEach(v => {nodeColors[v] = redColor});
    drawGraphImage(_global.graph, '#reductionRuleRedStarReductionMainGraphImage', nodeColors);
}

function renderReductionRuleRedComponentReductionMainGraphImage() {
    const node = document.getElementById('reductionRuleRedComponentReductionMainGraphImage');
    if(!node) return;
    const nodeColors = {};
    d3.range(_global.graph.n).forEach(v => {nodeColors[v] = '#aaa'});
    _global.rule.red_component1.forEach(v => {nodeColors[v] = redColor});
    _global.rule.red_component2.forEach(v => {nodeColors[v] = 'orange'});
    nodeColors[_global.rule.v] = 'purple'
    drawGraphImage(_global.graph, '#reductionRuleRedComponentReductionMainGraphImage', nodeColors);
}

function reductionRuleDetail() {
  renderReductionRuleRedStarReductionMainGraphImage();
  renderReductionRuleRedComponentReductionMainGraphImage();
}

function hookBranchingRuleTypeSelect() {
  const select = document.getElementById('branchingRuleTypeSelect');
  if(!select) return;
  select.addEventListener('change', event => {
    const value = event.target.value;
    document.querySelectorAll("tr[data-branching-rule-type]").forEach(elem => elem.classList.remove('d-none'))
    if(value!=='-') {
      document.querySelectorAll("tr[data-branching-rule-type]:not([data-branching-rule-type='"+value+"'])").forEach(elem => elem.classList.add('d-none'))
    }
  })
}

function init() {
  hookBranchingRuleTypeSelect();
  drawGraphWithRedVertices(_global.graph, '#mainGraphImage');
  reductionRuleDetail();

  for(const modal of document.querySelectorAll(".modal")) {
    modal.addEventListener('show.bs.modal', event => {
      console.log(modal.dataset)
      if(modal.dataset.subsetIdx) {
        renderModal({
          type: 'subset',
          render_id: 'subset_'+modal.dataset.subsetIdx,
          subset_idx: +modal.dataset.subsetIdx,
          subset: _global.rule.subsets[+modal.dataset.subsetIdx],
        })
      }
      if(modal.dataset.expansionIdx) {
        renderModal({
          type: 'expansion',
          render_id: 'expansion_'+modal.dataset.expansionIdx,
          expansion_idx: +modal.dataset.expansionIdx,
          expansion: _global.expansions[+modal.dataset.expansionIdx],
        })
      }
    })
  }
}



init();
