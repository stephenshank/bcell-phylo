import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';

const d3 = require('d3');
require("phylotree");


class Trees extends Component {
  constructor(props) {
    super(props);
    this.state = { newick: '' };
  }
  componentDidMount() {
    d3.text('/data/full_size-30.new', (err, data) => {
      const fills = { 6: 'pink', 7: 'lightblue', 9: 'plum', 11:'LemonChiffon '},
        strokes = { 6: 'red', 7: 'blue', 9: 'purple', 11: 'GoldenRod'};
      var tree = d3.layout.phylotree()
        .svg(d3.select("#tree_display"))
        .options({
          'left-right-spacing': 'fixed-step',
          'top-bottom-spacing': 'fixed-step'
        }).style_nodes(function(container, node){
          if(node.name[0]=='s') {
            const r = Math.floor(Math.sqrt(+node.name.split('_')[2].split('-')[1])),
              time = +node.name.split('_')[1].split('-')[1];
            container.append('circle')
              .attr('cx', r)
              .attr('cy', 0)
              .attr('r', r)
              .style('fill', fills[time])
              .style('stroke', strokes[time]);
            container.select('text')
              .attr('dx', 2*r+1)
              .attr('fill', strokes[time]);
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
      tree.layout();
    });
  }
  render(){
    return (<Row>
      <Col xs={12}>
        <svg id="tree_display" width={1200} height={1500}></svg>          
        <div style={{position: 'fixed', top: '70px', marginRight: '5px', marginTop: '5px', left: '1000px'}}>
          <svg>
            <rect width={20} height={20} fill={"red"} transform={'translate(0, 20)'} />
            <text x={25} y={35}>Time 6</text>
            <rect width={20} height={20} fill={"blue"} transform={'translate(0, 45)'} />
            <text x={25} y={60}>Time 7</text>
            <rect width={20} height={20} fill={"purple"} transform={'translate(0, 70)'} />
            <text x={25} y={85}>Time 9</text>
            <rect width={20} height={20} fill={"GoldenRod"} transform={'translate(0, 95)'} />
            <text x={25} y={110}>Time 11</text>
          </svg>
        </div>
      </Col>
    </Row>);
  }
}

module.exports = Trees;
