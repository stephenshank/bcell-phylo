import React, { Component } from 'react';
import { Row, Col } from 'react-bootstrap';
import Alignment from 'alignment.js';

const d3 = require('d3');


class Alignments extends Component {
  constructor(props){
    super(props);
    this.state = { fasta: '' };
  }
  componentDidMount() {
    d3.text('/data/full_size-30.fasta', (err, data) => {
      this.setState({ fasta: data });
    });
  }
  render() {
    return (<Row>
      <Col xs={12}>
        <Alignment fasta={this.state.fasta} width={1200} height={1500} />
      </Col>
    </Row>);
  }
}

module.exports = Alignments;
