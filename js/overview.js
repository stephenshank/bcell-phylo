import React, { Component } from 'react';
import ReactTable from "react-table";
import { withRouter } from "react-router-dom";
const d3 = require("d3");
import 'react-table/react-table.css'

const PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"];

const ExpansionTable = withRouter(function(props){
  const { history } = props;
  return (<ReactTable
    data={props.data}
    defaultSortDesc={true}
    className='-highlight'
    columns={[
      { Header: "Patient ID", accessor: 'patient_id', headerStyle: {fontWeight: 'bold'}},
      { Header: "V-gene allele", accessor: 'vgene_allele', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 1 size", accessor: 'visit1', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 2 size", accessor: 'visit2', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 3 size", accessor: 'visit3', headerStyle: {fontWeight: 'bold'}}
    ]}
    getTrProps={(state, rowInfo, column) => {
      const { patient_id, vgene_allele } = rowInfo.original,
        route = `/${patient_id}&${vgene_allele}`;
      return {
        onClick: ()=>history.push(route)
      }
    }}
  />);
});

class Overview extends Component {
  constructor(props) {
    super(props);
    this.state = {
      data: null,
      patient_id: 'all' 
    };
  }
  componentDidMount() {
    d3.csv('data/cluster_information.csv', (error, data) => {
      data.forEach(row=>{
        row.visit1 = +row.visit1;
        row.visit2 = +row.visit2;
        row.visit3 = +row.visit3;
      });
      this.setState({data: data});
    });
  }
  render() {
    if (!this.state.data) return <div/>;
    const data = this.state.patient_id == 'all' ?
      this.state.data :
      this.state.data.filter(datum => datum.patient_id == this.state.patient_id);
    return (<div style={{display: "flex", justifyContent: "center"}}>
      <div>
        <b> FILTERS - </b>
        <label>
          Patient ID: 
          <select value={this.state.patient_id} onChange={e=>this.setState({patient_id: e.target.value})}>
            <option value='all'>All</option>
            {PATIENT_IDS.map(patient_id => <option value={patient_id}>{patient_id}</option>)}
          </select>
        </label>
        <ExpansionTable data={data} /> 
      </div>
    </div>);
  }
}

module.exports = Overview;

