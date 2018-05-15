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
    this.loadJSON(7);
  }
  loadJSON(json_id) {
    d3.json(`/data/input/${json_id}_clone.json`, (err, data) => {
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
            defaultValue={7}
            onChange={value=>this.loadJSON(value)}
          >
            {[7, 8, 9, 10, 11, 12].map(d=>{
              return <ToggleButton value={d}>{d}</ToggleButton>
            })}
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
