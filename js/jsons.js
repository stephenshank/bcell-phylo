import React, { Component } from 'react';
import { ButtonToolbar, ToggleButtonGroup, ToggleButton, Row, Col } from 'react-bootstrap';
import ReactJson from 'react-json-view'

const d3 = require('d3');


class JSONs extends Component {
  constructor(props) {
    super(props);
    this.state = { 
      json_content: {}
    };
  }
  componentDidMount(){
    this.loadJSON(6);
  }
  loadJSON(json_id) {
    d3.json(`/data/${json_id}_clone.json`, (err, data) => {
      this.setState({ json_content: data });
    });
  }
  render() {
    return (<Row>
      <Col xs={12}>
        <ButtonToolbar>
          <ToggleButtonGroup
            type="radio"
            name="options"
            defaultValue={6}
            onChange={value=>this.loadJSON(value)}
          >
            <ToggleButton value={6}>6</ToggleButton>
            <ToggleButton value={7}>7</ToggleButton>
            <ToggleButton value={9}>9</ToggleButton>
            <ToggleButton value={11}>11</ToggleButton>
          </ToggleButtonGroup>
        </ButtonToolbar>
      </Col>
      <Col xs={12}>
        <ReactJson
          src={this.state.json_content}
          groupArraysAfterLength={5}
        />
      </Col>
    </Row>);
  }
}

module.exports = JSONs;
