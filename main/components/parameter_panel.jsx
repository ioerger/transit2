/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/events.js"
import { Column, Row, Input, Button, Code, askForFiles } from "../elements.jsx"
import { EasyFilePicker } from "./file_picker.jsx"
import { data, events } from "../data.jsx"

export const ParameterPanel = ({ children, style, }) => {
    return <Column width="100%" flex-grow="1">
        <span class="custom-header">
            Analysis Parameters
        </span>
        <Column
            padding="1.2rem 1rem"
            horizontalAlignment={`space-between`}
            verticalAlignment={"center"}
            background="whitesmoke"
            border="lightgray solid 1px"
            min-width="fit-content"
            border-radius="0.7rem"
            width="100%"
            flex-grow="1"
            gap="2rem"
            >
                <Row width="100%" horizontalAlignment="space-between">
                    Pseudocount <Input type="number" value="5" />
                </Row>
                
                <Row width="100%" horizontalAlignment="center">
                    <Button onClick={()=>alert("Running ...")} background="rgb(202, 104, 104)" color="white">
                        Run
                    </Button>
                </Row>
        </Column>
    </Column>
}

