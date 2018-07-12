import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';

const d3 = require('d3');
const _ = require('underscore');

import { BaseAlignment, fastaParser } from 'alignment.js';


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
    const { fasta, json } = nextProps;
    if(json && fasta) {
      this.displayTree(fasta, json);
    }
  }
  displayTree(fasta, json) {
    d3.select('#tree_display').html('');
    const size = 10*json.number_of_sequences;
    var tree = d3.layout.phylotree()
      .svg(d3.select("#tree_display"))
      .options({
        'left-right-spacing': 'fit-to-size',
        'top-bottom-spacing': 'fit-to-size',
        'selectable': false,
        'show-scale': false
      })
      .size([size, size]);
    tree(json.newick);
    console.log('size', tree.get_nodes().filter(d3.layout.phylotree.is_leafnode).length, json.number_of_sequences);
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
    this.sequence_data = fastaParser(fasta)
    this.sequence_data.sort((a,b)=>{
      const a_index = ordered_leaf_names.indexOf(a.header),
        b_index = ordered_leaf_names.indexOf(b.header);
      return a_index - b_index;
    }); 
  }
  render(){
    const style = {
      display: 'grid',
      gridTemplateColumns: '200px 500px 200px 300px',
      gridTemplateRows: '1000px'
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
        <div style={style}>

          <div>
            <svg height={500}>
              {legend.map((d,i) => {
                return [<rect width={20} height={20} fill={d.color} transform={`translate(0,${20+25*i})`} />,
                  <text x={25} y={35+25*i}>{d.text}</text>];
              }) } 
            </svg>
          </div>

          <div style={{width: 500, height: 1000, overflowX: 'scroll', overflowY: 'scroll'}}>
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

module.exports = Trees;
