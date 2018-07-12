import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import Trees from './trees.js';
import JSONs from './jsons.js';

require('./main.css');

class App extends Component {
  constructor(props){
    super(props);
    this.patients = [ '77612', '67029' ];
    this.genes = [1,2,3,4,5,6,7].map(i=>'V'+i);
    this.state = { 
      patient: null,
      gene: null
    };
  }
  loadData(patient, gene) {
    const json_path = `/data/out/${patient}/${gene}.json`;
    const fasta_path = `/data/out/${patient}/${gene}_aligned.fasta`;
    d3.text(fasta_path, (err, fasta_data) => {
      d3.json(json_path, (err, json_data) => {
        this.setState({
          patient: patient,
          gene: gene,
          fasta: fasta_data,
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
            <NavDropdown title='Patient'>
              {this.patients.map(patient_id => {
                const eventKey = { type: 'patient', value: patient_id };
                return (<MenuItem eventKey={eventKey} active={this.state.patient == patient_id}>
                  {patient_id}
                </MenuItem>);
                }) }
            </NavDropdown>
            <NavDropdown title='Gene'>
              {this.genes.map(gene => {
                const eventKey = { type: 'gene', value: gene };
                return (<MenuItem eventKey={eventKey} active={this.state.gene == gene}>
                  {gene}
                </MenuItem>);
              }) }
            </NavDropdown>
          </Nav>
        </Navbar>
        <Grid fluid>
          <Route exact path="/" render={ () => <Trees json={this.state.json} fasta={this.state.fasta} /> } />
        </Grid>
      </div>
    </Router>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)
