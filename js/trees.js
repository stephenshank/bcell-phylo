import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';

const d3 = require('d3');
const $ = require("jquery");
const _ = require('underscore');

import { BaseAlignment, fastaParser, ScrollBroadcaster } from 'alignment.js';


require("phylotree");


function get_name(node) {
  return node.name.split('_')[0];
}

function get_size(node) {
  return +node.name.split('_')[2].split('-')[1];
}

function get_time(node) {
  return +node.name.split('_')[1].split('-')[1];
}

function get_cdr3(node) {
  return node.name.split('_')[3];
}

class Trees extends Component {
  constructor(props) {
    super(props);
  }
  componentWillUpdate(nextProps) {
    const { json } = nextProps;
    if(json) {
      this.displayTree(json);
    }
  }
  componentDidUpdate() {
    const { scroll_broadcaster } = this;
    $("#tree-div").off("wheel");
    $("#tree-div").on("wheel", function(e) {
      scroll_broadcaster.handleWheel(e, 'freeze_x');
    });
  }
  displayTree(json) {
    d3.select('#tree_display').html('');
    const { number_of_sequences } = json;
    const size = this.props.site_size*number_of_sequences;
    var tree = d3.layout.phylotree()
      .svg(d3.select("#tree_display"))
      .options({
        'left-right-spacing': 'fit-to-size',
        'top-bottom-spacing': 'fit-to-size',
        'selectable': false,
        'show-scale': false,
        'align-tips': true
      })
      .size([size, size]);
    var parsed = d3.layout.newick_parser(json.newick);
    tree(parsed);
    tree.traverse_and_compute (function (node) {
      var d = 1;
      if (node.children && node.children.length) {
        d += d3.max (node.children, function (d) { return d["count_depth"];});
      }
      node["count_depth"] = d;
    });
    tree.resort_children (function (a,b) {
      return (a["count_depth"] - b["count_depth"]);
    });
    tree.node_circle_size(0).layout();

    const ordered_leaf_names = tree
      .get_nodes()
      .filter(d3.layout.phylotree.is_leafnode)
      .map(d=>d.name);
    this.sequence_data = fastaParser(json.fasta);
    const number_of_sites = this.sequence_data[0].seq.length;
    this.sequence_data.sort((a,b)=>{
      const a_index = ordered_leaf_names.indexOf(a.header),
        b_index = ordered_leaf_names.indexOf(b.header);
      return a_index - b_index;
    }); 

    var guide_tree = d3.layout.phylotree()
      .svg(d3.select('#guide-tree'))
      .options({
        'left-right-spacing': 'fit-to-size',
        'top-bottom-spacing': 'fit-to-size',
        'collapsible': false,
        'transitions': false,
        'show-scale': false,
        'brush': false,
        'selectable': false
      })
      .size([this.props.guide_size, this.props.guide_size])
      .node_circle_size(0);
    guide_tree(parsed).layout();

    var guide_scale = d3.scale.linear()
      .domain([0, size])
      .range([0, this.props.guide_size]);
    var rect = d3.select("#guide-tree")
      .append('rect')
        .attr('x', 12)
        .attr('y', 0)
        .attr('fill', 'red')
        .style('opacity', .8)
        .attr('width', guide_scale(this.props.tree_view_width))
        .attr('height', guide_scale(this.props.tree_view_height));

    this.scroll_broadcaster = new ScrollBroadcaster(
      { 
        width: number_of_sites*this.props.site_size,
        height: number_of_sequences*this.props.site_size 
      },
      {
        width: this.props.gridTemplateColumns[1], 
        height: this.props.gridTemplateRows[0]
      },
      { x_pixel: 0, y_pixel: 0 },
      [
        "alignmentjs-alignment"
      ]
    );

/*
    d3.select('#guide').on("click", function() {
      var coords = d3.mouse(this),
        new_x = x.invert(coords[0])-window.innerWidth/2,
        new_y = y.invert(coords[1])-window.innerHeight/2;
      window.scrollTo(new_x, new_y);
    });  
*/
  }
  render(){
    const style = {
      display: 'grid',
      gridTemplateColumns: this.props.gridTemplateColumns.join('px ')+'px',
      gridTemplateRows: this.props.gridTemplateRows.join('px ')+'px',
    },
    tree_style = {
      width: this.props.tree_view_width,
      height: this.props.tree_view_height,
      overflowX: 'scroll',
      overflowY: 'scroll'
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
        <div id='main_viz' style={style}>

          <div>
            <svg width={300} height={200}>
              {legend.map((d,i) => {
                return [<rect width={20} height={20} fill={d.color} transform={`translate(0,${20+25*i})`} />,
                  <text x={25} y={35+25*i}>{d.text}</text>];
              }) } 
            </svg>
            <svg width={this.props.guide_size} height={this.props.guide_size} id='guide-tree' />
          </div>

          <div id="tree-div" style={tree_style}>
            <svg id="tree_display" width={1200}></svg>
          </div>

          <div>

          </div>

          <div>
            <BaseAlignment
              width={300}
              height={1000}
              sequence_data={this.sequence_data}
              site_size={20}
            />
          </div>

        </div>
      </Col>
    </Row>);
  }
}

Trees.defaultProps = {
  'guide_size': 300,
  'tree_view_width': 500,
  'tree_view_height': 1000,
  'site_size': 20,
  'gridTemplateColumns': [300, 500, 0, 300],
  'gridTemplateRows': [1000]
}

module.exports = Trees;
