import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/string.js"

export const events = {
    annotationAdded: new Event(),
    wigDataAdded: new Event(),
}

export const data = {
    annotation: null,
    sampleFiles: null,
    conditions: null,
    panelInfo: null,
}