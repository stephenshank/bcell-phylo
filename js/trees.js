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
      var tree = d3.layout.phylotree()
        .svg(d3.select("#tree_display"));
      tree(data)
        .options({
          'left-right-spacing': 'fixed-step',
          'top-bottom-spacing': 'fixed-step'
        })
        .layout();
    });
  }
  render(){
    return (<Row>
      <Col xs={12}>
        <svg id="tree_display" width={1200} height={1500}></svg>          
      </Col>
    </Row>);
  }
}

module.exports = Trees;
