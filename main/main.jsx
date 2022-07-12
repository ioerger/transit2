/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "./elemental.js"
import { Column, Row, Input, askForFiles } from "./elements.jsx"

import { Annotation } from "./components/annotation.jsx"
import { WigLoader } from "./components/wig_loader.jsx"
import { SampleFileTable } from "./components/sample_file_table.jsx"
import { ConditionsTable } from "./components/conditions_table.jsx"
import { ParameterPanel } from "./components/parameter_panel.jsx"

// TODO:
    // dropdown menus
    // disabled checkboxes
    // status menu at the bottom

document.body.append(<Row height="100vh" width="100vw" overflow="hidden">
    <style>{`
        .custom-header {
            font-weight: 500;
            font-size: 14pt;
            color: gray;
            margin-bottom: 2px;
            margin-left: 2px;
        }
    `}</style>
    
    <Column name="FilePanel" width="27rem" background="white" padding="3rem" height="100%" max-height="100%">
        <Annotation />
        <Row height={`3rem`} />
        <Row flex-grow="1" min-height="20rem" height="100%">
            <WigLoader />
        </Row>
    </Column>
    
    <Column name="DataPanel" width="27rem" background="white" padding="3rem" height="100%" max-height="100%" flex-grow="1">
        <SampleFileTable />
        <Row height={`3rem`} />
        <ConditionsTable />
    </Column>
    
    <Column name="ParameterPanel" width="27rem" background="white" padding="3rem" height="100%" max-height="100%">
        <ParameterPanel />
    </Column>
</Row>)