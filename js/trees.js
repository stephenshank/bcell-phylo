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
          if(node.name[0]=='s') {
            const r = Math.floor(Math.sqrt(+node.name.split('_')[2].split('-')[1])),
              time = +node.name.split('_')[1].split('-')[1];
            container.append('circle')
              .attr('cx', r)
              .attr('cy', 0)
              .attr('r', r)
              .style('fill', fills[time])
              .style('stroke', 'black');
            container.select('text')
              .attr('dx', 2*r+1)
              .attr('fill', fills[time]);
          }
        }).node_span(function(node){
          if(d3.layout.phylotree.is_leafnode(node)) {
            return Math.floor(.5*(+node.name.split('_')[2].split('-')[1])**.5);
          }
        });
      tree(data);
      tree.traverse_and_compute (function (n) {
        var d = 1;
        if (n.children && n.children.length) {
          d += d3.max (n.children, function (d) { return d["count_depth"];});
        }
        n["count_depth"] = d;
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
                <text x={25} y={35+25*i}>Datapoint {d[0]}</text>];
            }) } 
          </svg>
        </div>
      </Col>
    </Row>);
  }
}

module.exports = Trees;
