/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "./elemental.js"
import { Column, Row, Input, askForFiles } from "./elements.jsx"

import { Annotation } from "./components/annotation.jsx"
import { WigLoader } from "./components/wig_loader.jsx"

document.body.append(<Row height="100vh">
    <Column name="GlobalPanel" width="70vw" background="white" padding="3rem" height="100%" max-height="100%" overflow="hidden">
        <style>{`
            .custom-header {
                font-weight: 500;
                font-size: 14pt;
                color: gray;
                margin-bottom: 2px;
                margin-left: 2px;
            }
        `}</style>
        <Annotation />
        <Row height={`3rem`} />
        <Row flex-grow="1" min-height="20rem" height="100%">
            <WigLoader />
        </Row>
        {/* <Column height="85vh">
            <Column padding="3rem">
                <button onClick={(event)=>console.log("file selected", event.target.files)} margin="2rem" >
                    Select File
                </button>
                <Column id="wigTable" background="gray" padding="3rem">
                    Table
                </Column>
            </Column>
            <Column padding="3rem">
                <button>
                    Load Files
                </button>
                <Column id="wigTable" background="gray" padding="3rem">
                    Table
                </Column>
            </Column>
        </Column> */}
    </Column>
    <Row name="LocalPanel" background="cornflowerblue" style="flex-grow: 1; align-items: center; justify-content: center;">
        panel
    </Row>
</Row>)