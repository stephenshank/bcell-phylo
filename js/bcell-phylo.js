import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';
const _ = require('underscore');

const d3 = require("d3");

const $ = require("jquery");

import {
  BaseAlignment,
  fastaParser,
  ScrollBroadcaster,
  SequenceAxis
} from 'alignment.js';


require("phylotree");


const colors = ['red', 'pink', 'blue', 'lightblue', 'purple', 'plum'];

function get_name(node, key) {
  key = key || 'name';
  return node[key].split('_')[0];
}

function get_size(node, key) {
  key = key || 'name';
  return +node[key].split('_')[2].split('-')[1];
}

function get_time(node, key) {
  key = key || 'name';
  return +node[key].split('_')[1].split('-')[1]-1;
}

function get_cdr3(node, key) {
  key = key || 'name';
  return node[key].split('_')[3];
}

function get_tag(node, key) {
  key = key || 'name';
  return node[key].split('_')[3];
}

function get_signature(node, key) {
  key = key || 'name';
  return node[key].split('_')[4];
}

function LegendItem(props) {
  key = key || 'name';
  const { text, color } = props;
}

function edgeStyler(element, edge) {
  if(!d3.layout.phylotree.is_leafnode(edge.target)) return null;
  const time = get_time(edge.target),
    stroke = colors[time];
  element.style("stroke", stroke);
}

class BCellPhylo extends Component {
  constructor(props) {
    super(props);

    this.column_sizes = [700, 700, 200, 200, 700];
    this.row_sizes = [20, 20, 700];
  }
  componentWillUpdate(nextProps) {

    const { site_size } = nextProps;
    this.sequence_data = fastaParser(nextProps.json.fasta)
      .map(record => {
        if (record.header.indexOf('Gene_V') > -1) return record;
        record.size = get_size(record, 'header');
        record.time = get_time(record, 'header');
        record.cdr3 = get_cdr3(record, 'header');
        record.old_header = record.header;
        //record.header = get_signature(record, 'header');
        return record;
      });
    const number_of_sequences = this.sequence_data.length;
    this.tree_size = number_of_sequences * site_size;
    this.main_tree = d3.layout
      .phylotree()
      .options({
        "left-right-spacing": "fit-to-size",
        "top-bottom-spacing": "fit-to-size",
        "show-scale": false,
        "align-tips": true,
        //"show-labels": false,
        selectable: false
      })
      .style_edges(edgeStyler)
      .size([this.tree_size, this.tree_size])
      .node_circle_size(0);
    this.parsed = d3.layout.newick_parser(nextProps.json.newick);
    this.main_tree(this.parsed);

    var i = 0;
    this.main_tree.traverse_and_compute(function(n) {
      var d = 1;
      if (!n.name) {
        n.name = "Node" + i++;
      }
      if (n.children && n.children.length) {
        d += d3.max(n.children, function(d) {
          return d["count_depth"];
        });
      }
      n["count_depth"] = d;
    });

    this.main_tree.resort_children(function(a, b) {
      return a["count_depth"] - b["count_depth"];
    }, true);

    const ordered_leaf_names = this.main_tree
      .get_nodes(true)
      .filter(d3.layout.phylotree.is_leafnode)
      .map(d => d.name);

    this.sequence_data.sort((a, b) => {
      if(a.header.indexOf('Gene_V') > -1) return -1;
      if(b.header.indexOf('Gene_V') > -1) return 1;
    
      const a_index = ordered_leaf_names.indexOf(a.header),
        b_index = ordered_leaf_names.indexOf(b.header);
      return a_index - b_index;
    });
  }
  componentDidUpdate() {
    this.main_tree.svg(d3.select("#alignmentjs-largeTreeAlignment")).layout();

    const guide_height = this.row_sizes[2],
      guide_width = this.column_sizes[0];

    this.guide_tree = d3.layout
      .phylotree()
      .svg(d3.select("#alignmentjs-guideTree"))
      .options({
        "left-right-spacing": "fit-to-size",
        // fit to given size top-to-bottom
        "top-bottom-spacing": "fit-to-size",
        // fit to given size left-to-right
        collapsible: false,
        // turn off the menu on internal nodes
        transitions: false,
        // turn off d3 animations
        "show-scale": false,
        // disable brush
        brush: false,
        // disable selections on this tree
        selectable: false
      })
      .style_edges(edgeStyler)
      .size([guide_height, guide_width])
      .node_circle_size(0);

    this.guide_tree(this.parsed).layout();

    this.guide_x_scale = d3.scale
      .linear()
      .domain([0, this.tree_size])
      .range([0, guide_width]);
    this.guide_y_scale = d3.scale
      .linear()
      .domain([0, this.tree_size])
      .range([0, guide_height]);
    this.rect = d3
      .select("#alignmentjs-guideTree")
      .append("rect")
      .attr("x", 0)
      .attr("y", 0)
      .attr("id", "guide-rect")
      .style("opacity", 0.8)
      .style("stroke-width", "2px")
      .style("stroke", "GoldenRod")
      .style("fill", "yellow")
      .attr("width", this.guide_x_scale(guide_width))
      .attr("height", this.guide_y_scale(guide_height));
     
    const alignment_axis_width = this.column_sizes[4],
      alignment_axis_height  = this.row_sizes[0];

    d3.select("#alignmentjs-axis-div")
      .style("width", alignment_axis_width + "px")
      .style("height", alignment_axis_height + "px");

    const number_of_sites = this.sequence_data[0].seq.length,
      number_of_sequences = this.sequence_data.length,
      { site_size } = this.props,
      alignment_width = site_size * number_of_sites,
      alignment_height = site_size * number_of_sequences;

    var alignment_axis_scale = d3.scale.linear()
      .domain([1, number_of_sites])
      .range([site_size / 2, alignment_width - site_size / 2]);

    var alignment_axis_svg = d3.select("#alignmentjs-alignment-axis");
    alignment_axis_svg.html("");
    alignment_axis_svg.attr("width", alignment_width)
      .attr("height", alignment_axis_height);

    var alignment_axis = d3.svg.axis()
      .orient("top")
      .scale(alignment_axis_scale)
      .tickValues(d3.range(1, number_of_sites, 2));

    alignment_axis_svg
      .append("g")
      .attr("class", "axis")
      .attr("transform", `translate(0, ${alignment_axis_height - 1})`)
      .call(alignment_axis);

    const bar_width = this.column_sizes[3],
      bar_height = this.row_sizes[1];

    var bar_axis_svg = d3.select("#alignmentjs-bar-axis");
    bar_axis_svg.html("");
    bar_axis_svg.attr("width", bar_width)
      .attr("height", bar_height);

    const { sequence_data } = this;
    console.log(sequence_data);
    var bar_scale = d3.scale.linear()
      .domain([0, d3.max(sequence_data.slice(1).map(d=>d.size))])
      .range([0, bar_width]);

    bar_axis_svg.append("g")
      .attr("class", "axis")
      .attr("transform", `translate(0, ${alignment_axis_height - 1})`)
      .call(bar_scale);

    var bar_svg = d3.select("#alignmentjs-bar");
    bar_svg.attr("width", this.column_sizes[3])
      .attr("height", alignment_height);

    //debugger;
    bar_svg.selectAll('rect')
      .data(sequence_data.slice(1))
      .enter()
        .append('rect')
        .attr('x', 0)
        .attr('y', function(d,i) { return i*site_size; })
        .attr('width', function(d) { return bar_scale(d.size); })
        .attr('height', site_size)
        .attr('fill', function(d) { return colors[d.time]; });
      
  }
  render(){
    if (!this.props.json) {
      return <div />;
    }
    const container_style = {
      display: "grid",
      gridTemplateColumns: this.column_sizes.join("px ") + "px",
      gridTemplateRows: this.row_sizes.join("px ") + "px"
    },
    legend = [
      { text: 'Visit 1, replicate 1', color: 'red' },
      { text: 'Visit 1, replicate 2', color: 'pink' },
      { text: 'Visit 2, replicate 1', color: 'blue' },
      { text: 'Visit 2, replicate 2', color: 'lightblue' },
      { text: 'Visit 3, replicate 1', color: 'purple' },
      { text: 'Visit 3, replicate 2', color: 'plum' }
    ];
    return (<Row>
      <Col xs={12}>
        <div id='main_viz' style={container_style}>

          <div style={{gridArea: "1 / 1 / 3 / 4"}}>
            <svg width={1400} height={40}>
              {legend.map((d,i) => {
                return (<g key={d.text} transform={`translate(${20+175*i}, 0)`} >
                  <rect width={20} height={20} fill={d.color} />,
                  <text x={25} y={15}>{d.text}</text>
                </g>);
              }) } 
            </svg>
            <svg width={this.props.guide_size} height={this.props.guide_size} id='guide-tree' />
          </div>

        <div></div>

        <div className="tree-scale-bar" style={{overflowX: "scroll"}}>
          <svg id="alignmentjs-alignment-axis" />
        </div>

        <div id="alignmentjs-bar-axis-div">
          <svg id="alignmentjs-bar-axis" />
        </div>

        <BaseAlignment
          width={this.column_sizes[4]}
          height={this.row_sizes[1]}
          sequence_data={[this.sequence_data[0]]}
          id='germline'
        />

        <div id="alignmentjs-guideTree-div">
          <svg id="alignmentjs-guideTree" />
        </div>

        <div
          id="alignmentjs-largeTreeAlignment-div"
          style={{ overflowX: "scroll", overflowY: "scroll" }}
        >
          <svg id="alignmentjs-largeTreeAlignment" />
        </div>

        <SequenceAxis
          width={this.column_sizes[2]}
          height={this.row_sizes[2]}
          sequence_data={this.sequence_data.slice(1)}
        />

        <div id="alignmentjs-bar-div" style={{overflowY: "scroll"}}>
          <svg id="alignmentjs-bar" />
        </div>

        <BaseAlignment
          width={this.column_sizes[4]}
          height={this.row_sizes[2]}
          sequence_data={this.sequence_data.slice(1)}
        />

        </div>
      </Col>
    </Row>);
  }
}

BCellPhylo.defaultProps = {
  'site_size': 20,
  'gridTemplateColumns': [300, 500, 0, 300],
  'gridTemplateRows': [1000]
}

module.exports = BCellPhylo;
