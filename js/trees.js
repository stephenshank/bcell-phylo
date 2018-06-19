import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';

const d3 = require('d3');
const _ = require('underscore');

require("phylotree");


const colors = {
  'no_replicates': {
      7: 'red',
      9: 'blue',
      11: 'purple'
    },
  'replicates': {
      7: 'red',
      8: 'pink',
      9: 'blue',
      10: 'lightblue',
      11: 'purple',
      12: 'plum'
    }
  },
  labels = {
    7: 'Visit 1, Tube 2',
    8: 'Visit 1, Tube 4',
    9: 'Visit 2, Tube 2',
    10: 'Visit 2, Tube 4',
    11: 'Visit 3, Tube 2',
    12: 'Visit 3, Tube 4'
  };

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
    this.state = { newick: '' };
  }
  componentDidUpdate() {
    this.displayTree();
  }
  componentDidMount() {
    this.displayTree();
  }
  displayTree(dataset) {
    const { active } = this.props;
    d3.select('#tree_display').html('');
    d3.text(`/data/${active}/full_size-30.new`, (err, data) => {
      const fills = colors[active];
      var tree = d3.layout.phylotree()
        .svg(d3.select("#tree_display"))
        .options({
          'left-right-spacing': 'fixed-step',
          'top-bottom-spacing': 'fixed-step',
          'selectable': false
        }).style_nodes(function(container, node){
          if(d3.layout.phylotree.is_leafnode(node)) {
            const r = Math.floor(Math.sqrt(get_size(node)));
            const p = Math.random();
            const arc = d3.svg.arc()
              .outerRadius(r)
              .innerRadius(0),
              pie = d3.layout.pie()
                .value(function(d) { return d[1]; }),
              fan_g = container.selectAll(".fan")
                .data(pie(_.pairs(node.pie_data)))
                .enter()
                  .append("g")
                  .attr("class", "fan");
              fan_g.append("path")
                .attr("d", arc)
                .attr("fill", function(d) { return fills[d.data[0]]; })
                .attr("transform", `translate(${r},${0})`);
            container.selectAll("text")
              .attr("transform", `translate(${2*r},0)`);
          }
        }).node_span(function(node){
          if(d3.layout.phylotree.is_leafnode(node)) {
            return Math.floor(.5*(get_size(node))**.5);
          }
        });
      tree(data);
      tree.traverse_and_compute (function (node) {
        if(node.children && _.all(node.children.map(d3.layout.phylotree.is_leafnode)) ) {
          const cdr3_1 = get_cdr3(node.children[0]),
            cdr3_2 = get_cdr3(node.children[1]);
          if(cdr3_1 == cdr3_2) {
            const child1_size = get_size(node.children[0]),
              child2_size = get_size(node.children[1]),
              child1_name = get_name(node.children[0]),
              child2_name = get_name(node.children[1]),
              new_size = child1_size+child2_size;

            node.name = child1_name + "+" + child2_name + "_MERGED_size-"+new_size+"_"+cdr3_1;
            node.pie_data = _.mapObject(node.children[0].pie_data, function(val, key) {
              return val+node.children[1].pie_data[key];
            });
            node.children = undefined;
          }
        } else if (d3.layout.phylotree.is_leafnode(node)) {
          node.pie_data = _.mapObject(fills, function(val, key) {
            const time = get_time(node);
            return key == time ? get_size(node) : 0;
          });
        }
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
    });
  }
  render(){
    const { active } = this.props,
      fills = colors[active];
    return (<Row>
      <Col xs={12}>
        <svg id="tree_display" width={1200} height={1500}></svg>          
        <div style={{position: 'fixed', top: '100px', marginRight: '5px', marginTop: '5px', left: '5%'}}>
          <svg height={500}>
            {_.pairs(fills).map((d,i) => {
              return [<rect width={20} height={20} fill={d[1]} transform={`translate(0,${20+25*i})`} />,
                <text x={25} y={35+25*i}>{labels[d[0]]}</text>];
            }) } 
          </svg>
        </div>
      </Col>
    </Row>);
  }
}

module.exports = Trees;
