import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import BCellPhylo from './bcell-phylo.js';
import JSONs from './jsons.js';

require('./main.css');

class App extends Component {
  constructor(props){
    super(props);
    this.patients = [ "28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS" ];
    this.genes = [1,2,3,4,5,6].map(i=>'V'+i);
    this.state = { 
      patient: null,
      gene: null
    };
  }
  loadData(patient, gene) {
    this.setState({json: null}, function() {
      const json_path = `/data/out/${patient}/${gene}.json`;
      d3.json(json_path, (err, json_data) => {
        json_data.patient = patient;
        json_data.gene = gene;
        this.setState({
          patient: patient,
          gene: gene,
          json: json_data
        });
      });
    });
  }
  componentDidMount() {
    const patient = '77612';
    const gene = 'V3';
    this.loadData(patient, gene);
  }
  onSelect(key){
    const patient = key.type == 'patient' ? key.value : this.state.patient;
    const gene = key.type == 'gene' ? key.value : "V3";
    this.loadData(patient, gene);
  }
  render(){
    return(<Router>
      <div>
        <Navbar onSelect={key=>this.onSelect(key)} fixedTop>
          <Navbar.Header>
            <Navbar.Brand>
              VEG Immunology Visualization
            </Navbar.Brand>
          </Navbar.Header>
          <Nav>
            <NavDropdown title='Patient' id='patient'>
              {this.patients.map(patient_id => {
                const eventKey = { type: 'patient', value: patient_id };
                return (<MenuItem
                  eventKey={eventKey}
                  active={this.state.patient == patient_id}
                  key={patient_id}
                >
                  {patient_id}
                </MenuItem>);
                }) }
            </NavDropdown>
            <NavDropdown title='Gene' id='patient'>
              {this.genes.map(gene => {
                const eventKey = { type: 'gene', value: gene };
                return (<MenuItem
                  eventKey={eventKey}
                  active={this.state.gene == gene}
                  key={gene}
                >
                  {gene}
                </MenuItem>);
              }) }
            </NavDropdown>
          </Nav>
        </Navbar>
        <Grid fluid>
          <Route exact path="/" render={ () => <BCellPhylo json={this.state.json} /> } />
        </Grid>
      </div>
    </Router>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)
